"""helper_modules/gspread_retry.py against a fake gspread (gspread itself is not needed to run these tests)."""
import os, sys, types, unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import requests
from helper_modules import gspread_retry as GR


class FakeAPIError(Exception):
    def __init__(self, code):
        super().__init__(f'APIError {code}'); self.response = types.SimpleNamespace(status_code=code)


def fake_gspread(script):
    """A gspread whose Client.request plays `script`: a list of outcomes (exceptions are raised, anything else is returned)."""
    g = types.ModuleType('gspread'); c = types.ModuleType('gspread.client'); ex = types.ModuleType('gspread.exceptions')
    ex.APIError = FakeAPIError
    calls = []
    class Client:
        def request(self, method, endpoint, params=None, **kw):
            calls.append((method, endpoint))
            o = script.pop(0) if script else 'ok'
            if isinstance(o, Exception):
                raise o
            return o
    c.Client = Client; g.client = c; g.exceptions = ex
    return g, c, ex, Client, calls


class T(unittest.TestCase):
    def setUp(self):
        self.script = []
        g, c, ex, self.Client, self.calls = fake_gspread(self.script)
        self.saved = {k: sys.modules.get(k) for k in ('gspread', 'gspread.client', 'gspread.exceptions', 'gspread.http_client')}
        sys.modules.update({'gspread': g, 'gspread.client': c, 'gspread.exceptions': ex}); sys.modules.pop('gspread.http_client', None)
        self.sleeps, self.logs = [], []
        self.old_sleep, GR._sleep = GR._sleep, self.sleeps.append
        GR.install(log=self.logs.append)

    def tearDown(self):
        GR._sleep = self.old_sleep
        for k, v in self.saved.items():
            if v is None: sys.modules.pop(k, None)
            else: sys.modules[k] = v

    def call(self, method='get', endpoint='https://sheets.googleapis.com/v4/spreadsheets/ID'):
        return self.Client().request(method, endpoint, params={})

    def test_a_connection_reset_is_retried_with_growing_waits(self):
        self.script[:] = [requests.exceptions.ConnectionError('Connection reset by peer'), requests.exceptions.ConnectionError('again'), 'metadata']
        self.assertEqual(self.call(), 'metadata'); self.assertEqual(len(self.calls), 3)
        self.assertEqual(len(self.sleeps), 2); self.assertTrue(3.5 <= self.sleeps[0] <= 6.5 and 7.5 <= self.sleeps[1] <= 12.5, self.sleeps)         # about 5 s, then about 10 s (with jitter)
        self.assertTrue('[gspread-retry] GET spreadsheets/ID failed (ConnectionError' in self.logs[0] and 'retry 1/5' in self.logs[0])

    def test_it_gives_up_after_the_last_try_and_raises_the_original_error(self):
        self.script[:] = [requests.exceptions.ConnectionError('down')] * 10
        with self.assertRaises(requests.exceptions.ConnectionError):
            self.call()
        self.assertEqual((len(self.calls), len(self.sleeps)), (6, 5)); self.assertLessEqual(max(self.sleeps), 60 * 1.25)

    def test_http_errors_are_retried_only_when_they_are_transient(self):
        self.script[:] = [FakeAPIError(429), FakeAPIError(503), 'ok']
        self.assertEqual(self.call(), 'ok'); self.assertEqual(len(self.sleeps), 2)
        self.script[:] = [FakeAPIError(404)]
        with self.assertRaises(FakeAPIError):
            self.call()                                                                                                    # not found / permission: a retry cannot help
        self.assertEqual(len(self.sleeps), 2)

    def test_a_request_that_might_have_reached_google_is_not_repeated_unless_it_is_idempotent(self):
        append = 'https://sheets.googleapis.com/v4/spreadsheets/ID/values/Blumer!A1:append'
        structure = 'https://sheets.googleapis.com/v4/spreadsheets/ID:batchUpdate'
        values = 'https://sheets.googleapis.com/v4/spreadsheets/ID/values:batchUpdate'
        for ep in (append, structure):
            del self.calls[:]; self.script[:] = [requests.exceptions.ConnectionError('reset'), 'ok']
            with self.assertRaises(requests.exceptions.ConnectionError):
                self.call('post', ep)
            self.assertEqual(len(self.calls), 1, ep)                                                                       # a repeated append would add the row twice
        del self.calls[:]; self.script[:] = [requests.exceptions.ConnectionError('reset'), 'ok']
        self.assertEqual(self.call('post', values), 'ok'); self.assertEqual(len(self.calls), 2)                             # writing the same values again changes nothing
        del self.calls[:]; self.script[:] = [FakeAPIError(503), 'ok']
        self.assertEqual(self.call('post', append), 'ok'); self.assertEqual(len(self.calls), 2)                             # 503 / 429: Google did not process it
        self.script[:] = [FakeAPIError(500)]
        with self.assertRaises(FakeAPIError):
            self.call('post', append)                                                                                      # 500 on an append: unknown whether it was applied

    def test_installing_never_raises_even_when_gspread_looks_different(self):
        GR.uninstall(); del self.Client.request                                                   # a gspread without the method the retry wraps
        logs = []
        self.assertFalse(GR.install(log=logs.append)); self.assertTrue(any('not installed (AttributeError' in l for l in logs))
        self.assertFalse(GR.install(log=lambda *a: (_ for _ in ()).throw(OSError('log is broken'))))           # not even a failing log function can make it raise

    def test_installing_twice_does_nothing_and_uninstall_restores_the_original(self):
        self.assertFalse(GR.install(log=self.logs.append)); self.assertTrue(self.Client.request._gspread_retry)
        self.assertEqual(self.Client.request.__name__, 'request')
        self.assertTrue(GR.uninstall()); self.assertFalse(hasattr(self.Client.request, '_gspread_retry'))
        self.script[:] = [requests.exceptions.ConnectionError('x')]
        with self.assertRaises(requests.exceptions.ConnectionError):
            self.call()                                                                                                    # no retry once removed


if __name__ == '__main__':
    unittest.main()
