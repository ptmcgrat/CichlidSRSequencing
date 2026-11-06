"""Process ONE ENA run, in two stages that the scheduler runs in separate pools.  Replaces grabENA.py.

  python -m download_reads.worker <spec.json>          spec['stage'] = fetch | convert | all   (started by `download_reads run`)

  fetch    download both mates (Aspera, one session; https fallback), md5 and size verified; the FASTQs stay on disk for the convert stage
  convert  FastqToSam on the FASTQs (cleaning them only if FastqToSam rejects them) -> MarkIlluminaAdapters -> upload -> verify the upload on the cloud -> delete everything

Exit status (grabENA.py exited 0 after failed conversions, so failed runs were recorded as downloaded):
  0 = done, 10 = failed, a retry will not help, 11 = failed, may be transient, 12 = the server / network is pushing back (the scheduler lowers the download concurrency).
A JSON status file says which stage failed and why. Platforms other than ILLUMINA / ELEMENT no longer fail at upload; everything large is deleted in a `finally`, also after a failure."""
import json, os, shutil, subprocess, sys, time, traceback
from download_reads.sources import ena as E, sra as SF

EXIT = {'permanent': 10, 'transient': 11, 'throttle': 12, 'source_unavailable': 13}


class StageFailed(RuntimeError):
    def __init__(self, stage, msg, kind='permanent'):
        super().__init__(f'{stage}: {msg}')
        self.stage, self.kind = stage, kind


def _tail(b, n=1500):
    s = b.decode('utf-8', 'replace') if isinstance(b, (bytes, bytearray)) else str(b or '')
    return s[-n:]


def fastq_to_sam(fq1, fq2, out_bam, spec):
    return ['gatk', '--java-options', spec.get('java_options', '-Xmx8g'), 'FastqToSam', '--FASTQ', fq1, '--FASTQ2', fq2, '--READ_GROUP_NAME', spec['run'], '--TMP_DIR', spec['temp_dir'],
            '--OUTPUT', out_bam, '--SAMPLE_NAME', spec['sample'], '--LIBRARY_NAME', spec['library'], '--PLATFORM', spec['platform']]


def _write_status(spec, st):
    st['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    if spec.get('status_json'):
        tmpf = spec['status_json'] + '.tmp'
        with open(tmpf, 'w') as f:
            json.dump(st, f, indent=1)
        os.replace(tmpf, spec['status_json'])


def _fastq_paths(spec):
    if spec.get('fastq'):                                           # the SRA route writes <run>_1.fastq / <run>_2.fastq (uncompressed)
        return list(spec['fastq'])
    return [os.path.join(spec['temp_dir'], spec['sel'][k]['name']) for k in ('r1', 'r2')]


def fetch_stage(spec, download=E.download_run, log=print, sra=None):
    rid, tmp = spec['run'], spec['temp_dir']
    st = dict(run=rid, state='failed', stage='download', kind='transient', error='', started=time.strftime('%Y-%m-%dT%H:%M:%S'), download=None)
    try:
        os.makedirs(tmp, exist_ok=True)
        if spec.get('source') == 'sra':
            paths, stats = (sra or SF.fetch_sra)(spec, tmp, log=log)
            st.update(state='fetched', stage='fetched', kind='', download=stats, fastq=paths)
            return st
        _, _, stats = download(spec['sel'], tmp, ascp=spec.get('ascp'), key=spec.get('key'), rate=spec.get('ascp_rate', '300m'), both_in_one_session=spec.get('one_session', True),
                               timeout=spec.get('ascp_timeout'), stall_s=spec.get('stall_s', 90.0), log=log)
        stats.setdefault('source', 'ENA-' + stats.get('method', '?'))
        st.update(state='fetched', stage='fetched', kind='', download=stats)
    except E.EnaError as e:
        st.update(error=str(e)[-1500:], kind=e.kind)
    except OSError as e:
        st.update(error=f'{type(e).__name__}: {e}', kind='transient')
    except Exception:
        st.update(error=traceback.format_exc()[-1500:], kind='permanent')
    finally:
        if st['state'] != 'fetched' and st['kind'] == 'permanent':           # keep partial files after a transient failure (resume), not after a permanent one
            for p in _fastq_paths(spec):
                for f in (p, p + '.partial'):
                    try:
                        os.remove(f)
                    except OSError:
                        pass
        _write_status(spec, st)
    return st


def convert_stage(spec, fm, run=subprocess.run, sanitize=E.sanitize_pair, log=print):
    rid, tmp, out = spec['run'], spec['temp_dir'], spec['output_bam']
    metrics = out + '.metrics.txt'; temp_bam = f'{tmp}{rid}_temp.bam'
    st = dict(run=rid, state='failed', stage='FastqToSam', kind='permanent', error='', started=time.strftime('%Y-%m-%dT%H:%M:%S'), sanitized=None, timings={})
    clock = time.time; tm = st['timings']; t_ = clock()
    fq1, fq2 = _fastq_paths(spec)
    created = [temp_bam, out, metrics, fq1, fq2, fq1 + '.partial', fq2 + '.partial', os.path.join(tmp, rid + '.sra')]
    try:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        if not (os.path.exists(fq1) and os.path.exists(fq2)):
            raise StageFailed('FastqToSam', 'the downloaded FASTQ files are not on disk (fetch stage not run, or cleaned)', 'transient')
        r = run(fastq_to_sam(fq1, fq2, temp_bam, spec), capture_output=True)
        if r.returncode != 0:
            err = _tail(r.stderr)
            if not E.EMPTY_READ_ERRORS.search(err):
                raise StageFailed('FastqToSam', err)
            st['stage'] = 'sanitize'                                  # fallback: zero-length reads etc.
            f1, f2 = f'{tmp}fixed_{rid}_1.fastq', f'{tmp}fixed_{rid}_2.fastq'; created += [f1, f2]
            n, kept = sanitize(fq1, fq2, f1, f2)
            st['sanitized'] = dict(pairs=n, kept=kept); st['stage'] = 'FastqToSam (sanitized)'
            r = run(fastq_to_sam(f1, f2, temp_bam, spec), capture_output=True)
            if r.returncode != 0:
                raise StageFailed('FastqToSam (sanitized)', _tail(r.stderr))
        tm['ubam_s'] = round(clock() - t_, 1); t_ = clock()
        st['stage'] = 'MarkIlluminaAdapters'
        plat = spec['platform']
        if plat == 'ILLUMINA' or 'ELEMENT' in plat:
            r = run(['gatk', 'MarkIlluminaAdapters', '-I', temp_bam, '-O', out, '-M', metrics, '--TMP_DIR', tmp], capture_output=True)
            if r.returncode != 0:
                raise StageFailed('MarkIlluminaAdapters', _tail(r.stderr))
        else:
            shutil.move(temp_bam, out)
            open(metrics, 'w').write(f'adapters not marked: platform {plat}\n')
        tm['adapters_s'] = round(clock() - t_, 1); t_ = clock()
        st['stage'] = 'upload'
        try:
            fm.uploadData(out); fm.uploadData(metrics)                  # FileManager.uploadData runs `rclone check` with check=True: a mismatch or a network problem raises
        except Exception as e:
            raise StageFailed('upload', f'{type(e).__name__}: {str(e)[-300:]}', 'transient')
        names = set(fm.returnCloudFiles(os.path.dirname(out) + '/'))
        missing = [os.path.basename(p) for p in (out, metrics) if os.path.basename(p) not in names]
        if missing:
            raise StageFailed('upload', f'{missing} not found on the cloud after upload', 'transient')
        tm['upload_s'] = round(clock() - t_, 1)
        st.update(state='done', stage='done', kind='')
    except StageFailed as e:
        st.update(stage=e.stage, error=str(e)[-1500:], kind=e.kind)
    except E.EnaError as e:
        st.update(error=str(e)[-1500:], kind=e.kind)
    except OSError as e:
        st.update(error=f'{type(e).__name__}: {e}', kind='transient')
    except Exception:
        st.update(error=traceback.format_exc()[-1500:], kind='permanent')
    finally:
        for f in created:
            try:
                os.remove(f)
            except OSError:
                pass
        shutil.rmtree(os.path.join(tmp, rid + '.fqtmp'), ignore_errors=True)
        _write_status(spec, st)
    return st


def process_run(spec, fm, run=subprocess.run, download=E.download_run, sanitize=E.sanitize_pair, log=print):
    """Both stages in one go (tests, and the case where one process is wanted)."""
    st = fetch_stage(dict(spec, status_json=None), download, log)
    if st['state'] != 'fetched':
        _write_status(spec, st)
        return st
    return convert_stage(spec, fm, run, sanitize, log)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit('usage: python -m download_reads.worker <spec.json>')
    spec = json.load(open(argv[0])); stage = spec.get('stage', 'all')
    fm = None
    if stage in ('convert', 'all'):
        from cichlid_common.file_manager import FileManager as FM
        fm = FM()
    st = fetch_stage(spec) if stage == 'fetch' else (convert_stage(spec, fm) if stage == 'convert' else process_run(spec, fm))
    print(f"{st['run']}: {st['state']} ({st['stage']}){' ' + st['error'][:300] if st['error'] else ''}")
    return 0 if st['state'] in ('done', 'fetched') else EXIT.get(st.get('kind'), 11)


if __name__ == '__main__':
    sys.exit(main())
