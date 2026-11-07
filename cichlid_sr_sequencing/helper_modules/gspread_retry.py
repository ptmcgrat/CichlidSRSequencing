"""Retry transient Google errors in EVERY gspread call (connection reset, timeouts, HTTP 429 / 500 / 502 / 503 / 504), by wrapping the one method all gspread calls go through (gspread's `request`).

Why: on utaka an alignment run died in `addAlignmentRow` with `Connection reset by peer` while opening the sheet. FileManager.setAlignmentDatabase retried its write three times, but the reads (readAlignmentDatabase,
readSampleDatabase, ...) and everything else did not, and one network hiccup killed the whole process.

    from helper_modules import gspread_retry; gspread_retry.install()        # file_manager.py does this at import once patched (unit_scripts/patch_file_manager.py)

Policy: up to 6 tries, waits of about 5, 10, 20, 40, 60 s (with jitter), a line printed for each retry. A request that is NOT idempotent is never repeated after a connection error, because the first one may have
reached Google (a repeated append would add a row twice): POSTs are retried after a connection error only for the value reads and writes that are idempotent (values:batchGet / batchUpdate / clear); any POST is
retried after HTTP 429 / 503 (Google did not process it). Python 3.7 compatible (utaka's environment)."""
import functools, random, re, time

TRIES, BASE, MAX_WAIT = 6, 5.0, 60.0
RETRY_CODES = (429, 500, 502, 503, 504)
NOT_PROCESSED = (429, 503)                                                       # the request was not processed: safe to repeat even when it is not idempotent
IDEMPOTENT_POST = re.compile(r'/values[^/]*:(batchGet|batchUpdate|clear|batchClear)$')
_sleep = time.sleep


def _target():
    """The class that owns `request`: HTTPClient in gspread >= 5.10, Client before."""
    try:
        from gspread import http_client
        if hasattr(http_client.HTTPClient, 'request'):
            return http_client.HTTPClient
    except ImportError:
        pass
    from gspread import client
    return client.Client


def install(tries=TRIES, base=BASE, max_wait=MAX_WAIT, log=print):
    """Wrap gspread's request method. Returns True if installed now, False if it already was."""
    import requests
    import gspread
    cls = _target()
    orig = cls.request
    if getattr(orig, '_gspread_retry', False):
        return False
    net = (requests.exceptions.ConnectionError, requests.exceptions.Timeout, requests.exceptions.ChunkedEncodingError)

    @functools.wraps(orig)
    def request(self, method, endpoint, *args, **kwargs):
        m = str(method).lower()
        safe_to_repeat = m != 'post' or bool(IDEMPOTENT_POST.search(str(endpoint)))
        err = None
        for attempt in range(1, tries + 1):
            try:
                return orig(self, method, endpoint, *args, **kwargs)
            except net as e:
                if not safe_to_repeat:
                    raise
                err = e
            except gspread.exceptions.APIError as e:
                code = getattr(getattr(e, 'response', None), 'status_code', 0)
                if code not in RETRY_CODES or (not safe_to_repeat and code not in NOT_PROCESSED):
                    raise
                err = e
            if attempt == tries:
                raise err
            wait = min(max_wait, base * 2 ** (attempt - 1)) * (0.75 + random.random() / 2)
            log('[gspread-retry] %s %s failed (%s: %s); retry %d/%d in %.0f s' % (m.upper(), str(endpoint).split('/v4/')[-1][:70], type(err).__name__, str(err)[:80], attempt, tries - 1, wait))
            _sleep(wait)
        raise err

    request._gspread_retry = True
    request._orig = orig
    cls.request = request
    return True


def uninstall():
    cls = _target()
    orig = getattr(cls.request, '_orig', None)
    if orig is not None:
        cls.request = orig
        return True
    return False
