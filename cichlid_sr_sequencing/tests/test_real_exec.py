"""The REAL subprocess path (Runner, bash pipe with pipefail, stdout redirection, process-group kill), with tiny fake gatk / minimap2 / samtools / rclone programs on the PATH.
What it proves that the fake-executor tests cannot: a crashed minimap2 inside the alignment pipe fails the sample (the old code only looked at the last process of the pipe),
and Ctrl-C style stops kill everything a command started."""
import os, shutil, stat, subprocess, sys, tempfile, threading, time, unittest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helper_modules import seq_pipeline as sp
from test_bam_stats import STATS, IDXSTATS, PICARD

PY = '#!/usr/bin/env python3\n'
TOOLS = {
    'gatk': PY + '''
import sys, os
a = sys.argv[1:]
def opt(k): return a[a.index(k) + 1]
if 'SamToFastq' in a:
    sys.stdout.write('@r1\\nACGT\\n+\\nIIII\\n')
elif 'MergeBamAlignment' in a:
    sys.stdin.read(); open(opt('-O'), 'w').write('bam')
elif 'MarkDuplicates' in a:
    out = opt('-O'); open(out, 'w').write('bam'); open(out[:-4] + '.bai', 'w').write('bai'); open(opt('-M'), 'w').write(%r)
elif 'MergeSamFiles' in a:
    open(opt('-O'), 'w').write('bam')
''' % PICARD,
    'minimap2': PY + '''
import sys, os, time
if '--version' in sys.argv: print('2.28-fake'); sys.exit(0)
sys.stdin.read()
if os.environ.get('FAKE_MM_PIDFILE'): open(os.environ['FAKE_MM_PIDFILE'], 'w').write(str(os.getpid()))
if os.environ.get('FAKE_MM_SLEEP'): time.sleep(float(os.environ['FAKE_MM_SLEEP']))
if os.environ.get('FAKE_MM_FAIL'): sys.exit(137)
print('@HD\\tVN:1.6')
''',
    'samtools': PY + '''
import sys
a = sys.argv[1:]
cmd = a[0]
if cmd == 'stats': sys.stdout.write(%r)
elif cmd == 'idxstats': sys.stdout.write(%r)
elif cmd == 'view' and '-c' in a: print(100000000)
elif cmd == 'view': open(a[a.index('-o') + 1], 'w').write('slice')
elif cmd == 'index':
    pos = [x for x in a[1:] if not x.startswith('-') and not x.isdigit()]
    open(pos[1] if len(pos) > 1 else pos[0] + '.bai', 'w').write('bai')
''' % (STATS, IDXSTATS),
    'rclone': PY + 'import sys\n',
    'fake_split': PY + '''
import sys, json
bam, prefix = sys.argv[2], sys.argv[3]
for t in ('unmapped', 'discordant', 'inversion', 'duplication', 'clipped', 'chimeric'): open(f'{prefix}.{t}.bam', 'w').write('bam')
json.dump({'reads_seen': 1}, open(prefix + '.split_summary.json', 'w'))
''',
}


class TestRealProcesses(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.bin = os.path.join(self.tmp, 'bin'); os.makedirs(self.bin)
        for name, body in TOOLS.items():
            path = os.path.join(self.bin, name)
            with open(path, 'w') as f:
                f.write(body)
            os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
        self.old_path = os.environ['PATH']
        os.environ['PATH'] = self.bin + os.pathsep + self.old_path

    def tearDown(self):
        os.environ['PATH'] = self.old_path
        for k in ('FAKE_MM_FAIL', 'FAKE_MM_SLEEP', 'FAKE_MM_PIDFILE'):
            os.environ.pop(k, None)
        shutil.rmtree(self.tmp, ignore_errors=True)

    def pipeline(self, sample='SMP_1', n_ubams=1, **kw):
        stop = threading.Event()
        cfg = sp.Config(max_cores=16, max_ram_gb=200, samples_in_flight=2, align_threads=4, prefetch=1, python=os.path.join(self.bin, 'fake_split'), nice=0, **kw)
        budget = sp.ResourceBudget(16, 200, stop)
        runner = sp.Runner(budget, nice=0, stop=stop)
        os.makedirs(f'{self.tmp}/reads', exist_ok=True)
        path = sp.SamplePaths(sample, [f'{self.tmp}/reads/{sample}_{i}.bam' for i in range(n_ubams)], f'{self.tmp}/bam/{sample}', f'{self.tmp}/temp/{sample}')
        rows = []

        class Proc:
            pid = 1
            def poll(s): return 0

        def start(job):
            with open(job.path, 'w') as f:
                f.write('bam')
            return Proc()

        pipe = sp.SeqPipeline(cfg, runner, '/ref.fna', '/ref.mmi', 962_000_000, start, lambda p: (['rclone', 'copy', p.bam_dir, 'x'], ['rclone', 'check', p.bam_dir, 'x']), rows.append,
                              genome_version='G', log=lambda l: None, check_fn=lambda p: (os.path.exists(p), ''), download_kwargs=dict(kill_fn=lambda proc: None, poll=0.002, admit_poll=0.01, backoff=0.0))
        return pipe, path, rows

    def test_the_whole_pipeline_with_real_processes(self):
        pipe, path, rows = self.pipeline()
        text, done, failed = pipe.run_batch([path])
        self.assertEqual((done, failed), (['SMP_1'], {}), text)
        self.assertEqual(rows[0]['PrimaryReads'], 100_000_000)
        self.assertAlmostEqual(rows[0]['DuplicateFraction'], 0.09)

    def test_a_minimap2_crash_inside_the_pipe_fails_the_sample(self):
        os.environ['FAKE_MM_FAIL'] = '1'
        pipe, path, rows = self.pipeline()
        text, done, failed = pipe.run_batch([path])
        self.assertEqual((done, list(failed)), ([], ['SMP_1']), text)
        self.assertIn('exit code 137', failed['SMP_1'])
        self.assertEqual(rows, [])
        self.assertFalse(os.path.exists(path.sorted_bam))

    def test_without_pipefail_the_same_crash_would_have_gone_unnoticed(self):
        """The old code looked only at the last process of the pipe. This is why: MergeBamAlignment (the last process) exits 0 after the aligner died."""
        os.environ['FAKE_MM_FAIL'] = '1'
        cmd = sp.align_command('/r/u.bam', '/g/ref.mmi', '/g/ref.fna', self.tmp + '/', self.tmp + '/o.bam', 4, 4, 8, [self.tmp + '/a', self.tmp + '/b', self.tmp + '/c'])
        self.assertEqual(subprocess.run(['bash', '-c', cmd[-1]]).returncode, 0, 'a plain pipe reports only the last process')
        self.assertNotEqual(subprocess.run(cmd).returncode, 0, 'with pipefail the crash is reported')

    def test_a_stop_request_kills_the_whole_process_group(self):
        pidfile = self.tmp + '/mm.pid'
        os.environ['FAKE_MM_PIDFILE'], os.environ['FAKE_MM_SLEEP'] = pidfile, '60'
        pipe, path, rows = self.pipeline()
        threading.Timer(1.5, pipe.stop.set).start()
        t0 = time.time()
        _, done, failed = pipe.run_batch([path])
        self.assertLess(time.time() - t0, 20, 'the stop request must not wait for the 60 s sleep')
        self.assertEqual(failed, {})
        self.assertTrue(os.path.exists(pidfile))
        with open(pidfile) as f:
            pid = int(f.read())

        def alive(pid):                      # a killed process that has not been reaped yet is a zombie ('Z'): dead for our purposes
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                return False
            try:
                with open(f'/proc/{pid}/stat') as f:
                    return f.read().rsplit(')', 1)[1].split()[0] != 'Z'
            except OSError:
                return False
        deadline = time.time() + 10
        while alive(pid) and time.time() < deadline:
            time.sleep(0.1)
        self.assertFalse(alive(pid), 'the fake minimap2 is still running after the stop request')

    def test_a_missing_program_is_a_clear_failure_not_a_crash(self):
        os.environ['PATH'] = self.tmp + '/nothing'
        pipe, path, rows = self.pipeline()
        pipe.runner.exec_fn = sp.real_exec
        text, done, failed = pipe.run_batch([path])
        self.assertEqual(list(failed), ['SMP_1'])


if __name__ == '__main__':
    unittest.main()
