"""bchrom's download sources: SRA (module 02's fetch_sra, faked here), ENA, and `auto` (SRA first, ENA as the fallback)."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yaml
from cichlid_common.common_light import config_path
from sd_genotypes.bchrom import ena, steps as S, workflow as W
from download_reads.sources.ena import EnaError

REPORT = ('run_accession\tfastq_ftp\tfastq_md5\tfastq_bytes\tbase_count\tread_count\tlibrary_layout\tsample_alias\n'
          'SRR9672201\tftp.sra.ebi.ac.uk/vol1/fastq/SRR967/001/SRR9672201/SRR9672201_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR967/001/SRR9672201/SRR9672201_2.fastq.gz\taa;bb\t5000;5000\t2000000\t109512218\tPAIRED\tx\n')


def cfg_with(source, **dl):
    cfg = yaml.safe_load(open(config_path('bchrom.yaml'))); cfg['download'].update(source=source, **dl); return cfg


def processor(cfg, sra=None, down=None, sleep=lambda s: None):
    work = tempfile.mkdtemp() + '/'
    rec = ena.parse_filereport(REPORT)['SRR9672201']
    return W.RunProcessor(dict(id='SRR9672201'), rec, work, cfg, dict(reference='/r', index='/i'), downloader=down or (lambda r, d: (os.makedirs(d, exist_ok=True), open(d + r['name'], 'wb').write(b'fq'), d + r['name'])[2]),
                          sra_fetcher=sra, sleep=sleep, log=lambda *a: None)


class FakeSra:
    def __init__(self, fail=()):
        self.fail, self.calls = list(fail), []
    def __call__(self, spec, dest, **kw):
        self.calls.append(dict(spec))
        if self.fail:
            kind = self.fail.pop(0)
            if kind:
                open(os.path.join(dest, 'partial.junk'), 'w').write('x')                  # a failed source leaves debris
                raise EnaError(f'fake {kind}', kind)
        ps = [os.path.join(dest, f'{spec["run"]}_{k}.fastq') for k in (1, 2)]
        [open(p, 'w').write('x') for p in ps]
        return ps, dict(source='SRA-aws', mb_per_s=100.0, quality_check='identical-to-ENA', count_check='exact')


class TestSources(unittest.TestCase):
    def test_sra_only_passes_ena_s_raw_read_count_and_the_first_file_url_and_finds_the_fastq_pair(self):
        sra = FakeSra(); rp = processor(cfg_with('sra'), sra)
        rp._fetch()
        s = sra.calls[0]
        self.assertEqual((s['ena_read_count'], s['paired'], s['fq_threads'], s['sra_connections']), (109512218, True, 8, 8))
        self.assertEqual(s['ena_url_r1'], 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR967/001/SRR9672201/SRR9672201_1.fastq.gz')
        self.assertEqual([os.path.basename(p) for p in rp._fastq_pair()], ['SRR9672201_1.fastq', 'SRR9672201_2.fastq'])

    def test_auto_falls_back_to_ena_for_a_run_sra_cannot_deliver_and_clears_the_debris(self):
        sra = FakeSra(fail=['source_unavailable']); rp = processor(cfg_with('auto'), sra)
        rp._fetch()
        self.assertEqual(len(sra.calls), 1)                                                                # source_unavailable is final for SRA: no retries
        self.assertEqual(sorted(os.listdir(rp.p.fastq)), ['SRR9672201_1.fastq.gz', 'SRR9672201_2.fastq.gz'])    # ENA files only: the SRA debris is gone

    def test_transient_sra_errors_are_retried_with_a_backoff(self):
        sra = FakeSra(fail=['transient', 'throttle']); waits = []
        rp = processor(cfg_with('sra', attempts=4, backoff_s=30), sra, sleep=waits.append)
        rp._fetch()
        self.assertEqual((len(sra.calls), waits), (3, [30.0, 60.0]))

    def test_sra_only_does_not_fall_back_and_says_why(self):
        rp = processor(cfg_with('sra'), FakeSra(fail=['source_unavailable']))
        with self.assertRaises(ena.DownloadError) as cm:
            rp._fetch()
        self.assertIn('source_unavailable', str(cm.exception))

    def test_all_sources_failing_reports_each(self):
        def down(r, d): raise ena.DownloadError('ENA is down')
        rp = processor(cfg_with('auto'), FakeSra(fail=['source_unavailable']), down=down)
        with self.assertRaises(ena.DownloadError) as cm:
            rp._fetch()
        self.assertTrue('sra:' in str(cm.exception) and 'ena: ENA is down' in str(cm.exception))

    def test_ena_only_is_unchanged(self):
        sra = FakeSra(); rp = processor(cfg_with('ena'), sra); rp._fetch()
        self.assertEqual(sra.calls, []); self.assertEqual(len(rp.files), 2)

    def test_disk_is_reserved_for_the_sra_route(self):
        self.assertEqual(processor(cfg_with('ena'), None)._disk_factor(), 5.0)
        self.assertEqual((processor(cfg_with('sra'), None)._disk_factor(), processor(cfg_with('auto'), None)._disk_factor()), (9.0, 9.0))

    def test_preflight_needs_the_sra_tools_only_when_they_are_used(self):
        have = lambda t: None if t in ('aria2c', 'fasterq-dump') else '/bin/' + t
        g = dict(reference='/r', index='/i')
        self.assertEqual(S.preflight(g, which=have, exists=lambda p: True, cfg=cfg_with('ena')), [])
        bad = S.preflight(g, which=have, exists=lambda p: True, cfg=cfg_with('auto'))
        self.assertTrue(any('aria2c' in b for b in bad) and any('fasterq-dump' in b for b in bad) and not any('vdb-validate' in b for b in bad))

    def test_a_bad_source_name_is_refused(self):
        with self.assertRaises(SystemExit):
            S.validate_config(cfg_with('prefetch'))


if __name__ == '__main__':
    unittest.main()
