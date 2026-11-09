"""The identity of the BAM in every panel JSON (`bam` = file name, `bam_size_bytes`): what is written, how a realigned BAM is told from the one measured, old (schema 1) files, the call tables, samples streamed
through retro (identity of the SOURCE bam), the HTTP size lookup against a real local server, and the command line exit codes."""
import functools, http.server, json, os, subprocess, sys, tempfile, threading, unittest, urllib.error
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import yaml
from cichlid_common.common_light import module_dirs, REPO_ROOT
from sd_genotypes.panel import identity as ID, measure as pm
from sd_genotypes.bchrom import measure as BM, score as BS
from sd_genotypes.lg5inv import call as LC
from sd_genotypes.te import call as TC
import test_bchrom_panel as TB, test_lg5inv_panel as TL, test_te_panel as TT, test_retro_parity as TR


def make_bam(d, name='S1.all.bam', size=1234):
    p = os.path.join(d, name); open(p, 'wb').write(b'x' * size); open(p + '.bai', 'wb').write(b'i'); return p


def bchrom_json(d, sample='S1', size=1234, name=None):
    cfg, _ = TB.make_def(); blocks, _, _ = BM.load_definition(cfg)
    bam = make_bam(d, name or f'{sample}.all.bam', size)
    res = pm.measure(bam, sample, modules=('bchrom',), run=TB.fake_run(cfg, TB.depth_fn(blocks, 10, 1.0)), cfgs={'bchrom': cfg})
    path = os.path.join(d, f'{sample}.panel.json'); pm.write_json(path, res); return cfg, bam, res, path


class TestWritten(unittest.TestCase):
    def test_the_json_carries_the_file_name_and_the_size_not_a_path(self):
        d = tempfile.mkdtemp(); cfg, bam, res, path = bchrom_json(d, size=1234)
        self.assertEqual((res['schema'], res['bam'], res['bam_size_bytes']), (2, 'S1.all.bam', 1234))
        self.assertNotIn('bam_size', res); self.assertNotIn(os.sep, res['bam'])
        self.assertEqual(json.load(open(path))['bam_size_bytes'], 1234)                       # and the file on disk says the same

    def test_an_identity_can_be_given_for_a_temporary_regional_bam(self):
        d = tempfile.mkdtemp(); cfg, _ = TB.make_def(); blocks, _, _ = BM.load_definition(cfg); bam = make_bam(d, 'S1.regions.bam', 10)
        res = pm.measure(bam, 'S1', modules=('bchrom',), run=TB.fake_run(cfg, TB.depth_fn(blocks, 10, 1.0)), cfgs={'bchrom': cfg}, identity=dict(bam='S1.all.bam', bam_size_bytes=9_000_000_000))
        self.assertEqual((res['bam'], res['bam_size_bytes']), ('S1.all.bam', 9_000_000_000))


class TestSameBam(unittest.TestCase):
    def test_a_realigned_bam_with_the_same_name_is_told_by_its_size(self):
        d = tempfile.mkdtemp(); cfg, bam, res, path = bchrom_json(d, size=1234)
        self.assertEqual(ID.same_bam(path, bam)[0], 'same'); self.assertEqual(ID.same_bam(res, bam)[0], 'same')   # a path or the parsed JSON
        make_bam(d, 'S1.all.bam', 1300)                                                        # realigned: same name, other size
        st, why = ID.same_bam(path, bam); self.assertEqual(st, 'different'); self.assertIn('1234', why); self.assertIn('1300', why); self.assertIn('realigned', why)
        other = make_bam(d, 'S1.new.bam', 1234)
        st, why = ID.same_bam(path, other); self.assertEqual(st, 'different'); self.assertIn('S1.all.bam', why)
        self.assertEqual(ID.same_bam(path, dict(bam='S1.all.bam', bam_size_bytes=1234))[0], 'same')              # an identity dict works too

    def test_old_schema_1_files_are_read_and_a_streamed_one_has_no_known_size(self):
        disk = dict(schema=1, sample='S', bam='/home/student/Temp/S.all.bam', bam_size=999, panel_version='v', measured='t', samtools='s', modules={})
        self.assertEqual(ID.normalize(disk), ('S.all.bam', 999))
        self.assertEqual(ID.same_bam(disk, dict(bam='S.all.bam', bam_size_bytes=999))[0], 'same'); self.assertEqual(ID.same_bam(disk, dict(bam='S.all.bam', bam_size_bytes=1000))[0], 'different')
        streamed = dict(disk, bam='/tmp/x/S.regions.bam', bam_size=12345, source='http://127.0.0.1:5572/S/S.all.bam')       # the regional BAM's path and size were stored: not the identity of the source
        self.assertEqual(ID.normalize(streamed), ('S.all.bam', None))
        st, why = ID.same_bam(streamed, dict(bam='S.all.bam', bam_size_bytes=12345)); self.assertEqual(st, 'unknown'); self.assertIn('does not record', why)

    def test_schemas_1_and_2_are_read_and_others_refused(self):
        d = tempfile.mkdtemp()
        for sch in (1, 2, 3):
            pm.write_json(f'{d}/s{sch}.panel.json', dict(schema=sch, sample=f's{sch}', bam='x.bam', bam_size=1, bam_size_bytes=1, panel_version='v', measured='t', samtools='s', modules={}))
        self.assertEqual(pm.read_panel_json(f'{d}/s1.panel.json')['schema'], 1); self.assertEqual(pm.read_panel_json(f'{d}/s2.panel.json')['schema'], 2)
        with self.assertRaises(ValueError):
            pm.read_panel_json(f'{d}/s3.panel.json')

    def test_collect_gives_every_row_the_identity_old_and_new(self):
        d = tempfile.mkdtemp(); cfg, bam, res, path = bchrom_json(d, 'NEW', size=500)
        old = dict(schema=1, sample='OLD', bam='/a/b/OLD.all.bam', bam_size=700, panel_version='v', measured='t', samtools='s', modules={'bchrom': dict(res['modules']['bchrom'])})
        pm.write_json(f'{d}/OLD.panel.json', old)
        rows, _ = pm.collect(d, 'bchrom'); by = {r['SampleID']: r for r in rows}
        self.assertEqual((by['NEW']['bam'], by['NEW']['bam_size_bytes'], by['OLD']['bam'], by['OLD']['bam_size_bytes']), ('NEW.all.bam', 500, 'OLD.all.bam', 700))
        self.assertNotIn('bam_size', by['NEW'])


class TestCallTables(unittest.TestCase):
    def test_every_call_table_carries_the_identity(self):
        ident = dict(bam='S1.all.bam', bam_size_bytes=4321)
        d = tempfile.mkdtemp(); cfg, bam, res, path = bchrom_json(d, size=4321)
        rows, _ = pm.collect(d, 'bchrom'); t = BS.call_table(rows, cfg)
        self.assertEqual((t.bam[0], int(t.bam_size_bytes[0])), ('S1.all.bam', 4321))
        tcfg, _ = TT.make_def(); out = TT.M.measure_bam('x.bam', tcfg, TT.fake_run(tcfg, {'chr7ins_c1': (4, 0), 'chr7ins_c2': (4, 0)}))
        t = TC.call_table([dict(SampleID='S1', **ident, **out)], tcfg)
        self.assertEqual((t.bam[0], int(t.bam_size_bytes[0])), ('S1.all.bam', 4321))
        raw = TL.raw(np.random.default_rng(1), p=0.5)
        t = LC.call_table([dict(SampleID='S1', **ident, **raw)], TL.CFG0)
        self.assertEqual((t.bam[0], int(t.bam_size_bytes[0])), ('S1.all.bam', 4321))


class TestRetro(unittest.TestCase):
    def test_a_streamed_sample_carries_the_identity_of_the_source_bam_not_of_the_temporary_one(self):
        url = 'http://fake/S1/S1.all.bam'
        res = TR.R.measure_sample('S1', url, ['gsdf'], tempfile.mkdtemp(), run=TR.fake_run(), head=TR.fake_head)
        self.assertEqual((res['bam'], res['bam_size_bytes'], res['source']), ('S1.all.bam', TR.fake_head(url), url))

    def test_a_missing_source_bam_is_no_bam_and_nothing_is_fetched(self):
        log = []
        with self.assertRaises(TR.R.FetchError) as e:
            TR.R.measure_sample('S3', 'http://fake/S3/S3.all.bam', ['gsdf'], tempfile.mkdtemp(), run=TR.fake_run(log), head=TR.fake_head)
        self.assertEqual((e.exception.kind, log), ('no_bam', []))

    def test_a_local_path_gives_its_own_name_and_size(self):
        p = make_bam(tempfile.mkdtemp(), 'L.all.bam', 77)
        self.assertEqual(ID.source_identity(p, head=lambda u: 1 / 0), dict(bam='L.all.bam', bam_size_bytes=77))
        self.assertEqual(ID.source_identity('http://h:1/a%20b/S%2D1.all.bam', head=lambda u: 5), dict(bam='S-1.all.bam', bam_size_bytes=5))


class _Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *a):
        pass


class _NoLength(_Quiet):
    """A server whose HEAD answers carry no Content-Length but which supports a one-byte range request."""
    def do_HEAD(self):
        self.send_response(200); self.send_header('Content-Type', 'application/octet-stream'); self.end_headers()

    def do_GET(self):
        n = os.path.getsize(self.translate_path(self.path))
        self.send_response(206); self.send_header('Content-Range', f'bytes 0-0/{n}'); self.send_header('Content-Length', '1'); self.end_headers(); self.wfile.write(b'x')


class TestHttpSize(unittest.TestCase):
    def serve(self, handler, d):
        srv = http.server.ThreadingHTTPServer(('127.0.0.1', 0), functools.partial(handler, directory=d))
        threading.Thread(target=srv.serve_forever, daemon=True).start(); self.addCleanup(srv.shutdown); self.addCleanup(srv.server_close)
        return f'http://127.0.0.1:{srv.server_address[1]}'

    def test_head_gives_the_size_and_a_missing_file_is_a_404(self):
        d = tempfile.mkdtemp(); os.makedirs(d + '/S1'); open(d + '/S1/S1.all.bam', 'wb').write(b'x' * 4096); base = self.serve(_Quiet, d)
        self.assertEqual(ID.http_size(base + '/S1/S1.all.bam'), 4096)
        self.assertEqual(ID.source_identity(base + '/S1/S1.all.bam'), dict(bam='S1.all.bam', bam_size_bytes=4096))
        with self.assertRaises(urllib.error.HTTPError) as e:
            ID.http_size(base + '/S2/S2.all.bam')
        self.assertEqual(e.exception.code, 404)

    def test_a_server_without_content_length_is_asked_for_one_byte(self):
        d = tempfile.mkdtemp(); open(d + '/f.bam', 'wb').write(b'x' * 999); base = self.serve(_NoLength, d)
        self.assertEqual(ID.http_size(base + '/f.bam'), 999)


class TestCli(unittest.TestCase):
    def test_exit_codes_current_stale_unknown(self):
        d = tempfile.mkdtemp(); cfg, bam, res, path = bchrom_json(d, size=1234)
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(module_dirs()))
        run = lambda j, b: subprocess.run([sys.executable, '-m', 'sd_genotypes.panel.identity', '--json', j, '--bam', b], capture_output=True, text=True, cwd=REPO_ROOT, env=env)
        r = run(path, bam); self.assertEqual(r.returncode, 0, r.stderr[-200:]); self.assertIn('CURRENT', r.stdout)
        make_bam(d, 'S1.all.bam', 2000)
        r = run(path, bam); self.assertEqual(r.returncode, 1); self.assertIn('STALE', r.stdout)
        legacy = d + '/legacy.panel.json'; pm.write_json(legacy, dict(schema=1, sample='L', bam='/x/y.bam', bam_size=5, source='http://h/y/y.all.bam', panel_version='v', measured='t', samtools='s', modules={}))
        r = run(legacy, bam); self.assertEqual(r.returncode, 2); self.assertIn('UNKNOWN', r.stdout)


if __name__ == '__main__':
    unittest.main()
