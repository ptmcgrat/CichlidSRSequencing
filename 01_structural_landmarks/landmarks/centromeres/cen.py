"""Centromere calling from satellite-repeat density. Pure functions (unit-tested). Cichlid centromeres are often poorly assembled, so a call is a satellite-density
peak with a confidence class, not a measured centromere; 'none' means not located."""
import gzip, re
import numpy as np
import pandas as pd


# ------------------------------------------------------------------ annotation parsing (streaming, keeps only wanted classes)
def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def _wanted(cls, classes):
    c = (cls or '').lower()
    return any(c.startswith(k.lower()) for k in classes)


def parse_annotation(path, classes):
    """RepeatMasker .out, GFF3, UCSC rmsk table or BED (tandem-repeat arrays; class assumed 'Satellite'). Returns DataFrame chrom,start,end,name,cls (0-based, half-open)."""
    rows = []
    with _open(path) as f:
        header = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith('#bin') or 'genoName' in line[:200]:
                header = line.lstrip('#').split(); continue
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue
            tab = line.rstrip('\n').split('\t')
            if header and len(tab) >= len(header) - 1 and 'genoName' in header:           # UCSC rmsk
                d = dict(zip(header, tab))
                cls = d.get('repClass', '') + '/' + d.get('repFamily', '')
                if _wanted(cls, classes):
                    rows.append((d['genoName'], int(d['genoStart']), int(d['genoEnd']), d.get('repName', ''), cls))
                continue
            if path.lower().replace('.gz', '').endswith(('.bed',)):
                if len(tab) >= 3:
                    rows.append((tab[0], int(tab[1]), int(tab[2]), tab[3] if len(tab) > 3 else 'tandem', 'Satellite'))
                continue
            if len(tab) == 9 and tab[3].isdigit() and tab[4].isdigit():                    # GFF3
                attrs = dict(kv.split('=', 1) for kv in tab[8].split(';') if '=' in kv)
                cls = next((attrs[k] for k in ('class_family', 'Class', 'class', 'repeat_class', 'family') if k in attrs), tab[2])
                name = (attrs.get('Target', '').split(' ')[0].replace('Motif:', '') or attrs.get('Name') or attrs.get('ID') or tab[2])
                if _wanted(cls, classes):
                    rows.append((tab[0], int(tab[3]) - 1, int(tab[4]), name, cls))
                continue
            sp = line.split()                                                               # RepeatMasker .out
            if len(sp) >= 11 and sp[0].isdigit() and sp[5].isdigit() and sp[6].isdigit():
                if _wanted(sp[10], classes):
                    rows.append((sp[4], int(sp[5]) - 1, int(sp[6]), sp[9], sp[10]))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'name', 'cls'])


# ------------------------------------------------------------------ FASTA lengths and N gaps (streaming)
def fasta_lengths_and_gaps(path):
    """dict name -> (length, [(start, end), ...] runs of N)."""
    out, name, pos, runs, cur = {}, None, 0, [], None
    with _open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    if cur is not None: runs.append(cur)
                    out[name] = (pos, runs)
                name, pos, runs, cur = line[1:].split()[0], 0, [], None
                continue
            s = line.strip()
            if 'N' in s or 'n' in s:
                for m in re.finditer(r'[Nn]+', s):
                    a, b = pos + m.start(), pos + m.end()
                    if cur is not None and cur[1] == a: cur = (cur[0], b)
                    else:
                        if cur is not None: runs.append(cur)
                        cur = (a, b)
            elif cur is not None and cur[1] < pos:
                runs.append(cur); cur = None
            pos += len(s)
        if name is not None:
            if cur is not None: runs.append(cur)
            out[name] = (pos, runs)
    return out


# ------------------------------------------------------------------ density and calling
def window_coverage(starts, ends, length, window):
    """Fraction of each window covered by the (merged) intervals."""
    n = int(np.ceil(length / window)); cov = np.zeros(n)
    if len(starts) == 0:
        return cov
    o = np.argsort(starts); s = np.asarray(starts)[o]; e = np.asarray(ends)[o]
    ms, me, cs, ce = [], [], s[0], e[0]
    for a, b in zip(s[1:], e[1:]):
        if a <= ce: ce = max(ce, b)
        else: ms.append(cs); me.append(ce); cs, ce = a, b
    ms.append(cs); me.append(ce)
    for a, b in zip(ms, me):
        w0, w1 = a // window, (b - 1) // window
        if w0 == w1: cov[w0] += b - a
        else:
            cov[w0] += (w0 + 1) * window - a; cov[w1] += b - w1 * window
            if w1 - w0 > 1: cov[w0 + 1:w1] += window
    return np.minimum(cov / window, 1.0)


def smooth(x, k):
    k = max(1, int(k))
    if k == 1: return np.asarray(x, float)
    return np.convolve(np.pad(np.asarray(x, float), k // 2, mode='edge'), np.ones(k) / k, mode='valid')[:len(x)]


def call_centromere(sat, window, smooth_bp, conf_cfg, gaps=None, floor=0.005, min_peak=0.02):
    """sat: per-window satellite fraction. Returns dict(start, end, mid, peak, ratio, dominance, confidence, gap_fraction)."""
    s = smooth(sat, smooth_bp // window); n = len(s)
    res = dict(start=np.nan, end=np.nan, mid=np.nan, peak=float(s.max()) if n else 0.0, ratio=np.nan, dominance=np.nan, confidence='none', gap_fraction=np.nan)
    if n == 0 or s.max() < min_peak:
        return res
    pk = int(np.argmax(s)); half = 0.5 * s[pk]
    lo = pk
    while lo > 0 and s[lo - 1] >= half: lo -= 1
    hi = pk
    while hi < n - 1 and s[hi + 1] >= half: hi += 1
    k = max(1, smooth_bp // window)
    outside = np.r_[s[:max(0, lo - k)], s[hi + 1 + k:]]
    second = outside.max() if len(outside) else 0.0
    ratio = s[pk] / max(np.median(s), floor); dom = s[pk] / max(second, floor)
    conf = 'low'
    for level in ('high', 'medium'):
        c = conf_cfg[level]
        if ratio >= c['min_ratio'] and dom >= c['min_dominance'] and s[pk] >= c.get('min_peak', 0.0):     # absolute floor: a faint bump is never 'high'
            conf = level; break
    res.update(start=int(lo * window), end=int((hi + 1) * window), mid=int((lo + hi + 1) * window // 2), ratio=float(ratio), dominance=float(dom), confidence=conf)
    if gaps is not None and len(gaps):
        res['gap_fraction'] = float(np.mean(gaps[lo:hi + 1]))
    return res


def top_families(sub, start, end, k=3):
    """Satellite families covering most bases inside the called interval: 'name:bp;...'."""
    if not len(sub):
        return ''
    inside = sub[(sub.end > start) & (sub.start < end)].copy()
    inside['bp'] = np.minimum(inside.end, end) - np.maximum(inside.start, start)
    t = inside.groupby('name').bp.sum().sort_values(ascending=False).head(k)
    return ';'.join(f'{n}:{int(b)}' for n, b in t.items())


# ------------------------------------------------------------------ array-level definition (family-specific, feature resolution)
def merge_arrays(starts, ends, max_gap=50_000):
    """Merge features of ONE satellite family into arrays: features closer than max_gap belong to the same array. Returns [(start, end, covered_bp, n_features)]."""
    if len(starts) == 0:
        return []
    o = np.argsort(starts); s = np.asarray(starts)[o]; e = np.asarray(ends)[o]
    arrays, cs, ce, cov, n = [], int(s[0]), int(e[0]), int(e[0] - s[0]), 1
    for a, b in zip(s[1:], e[1:]):
        a, b = int(a), int(b)
        if a - ce <= max_gap:
            cov += max(0, b - max(a, ce)); ce = max(ce, b); n += 1
        else:
            arrays.append((cs, ce, cov, n)); cs, ce, cov, n = a, b, b - a, 1
    arrays.append((cs, ce, cov, n))
    return arrays


def call_by_arrays(arrays, length, min_end_distance=500_000, min_array=100_000, ambiguity_ratio=0.5):
    """The centromere = the LARGEST array (by covered bases) that lies more than min_end_distance from BOTH chromosome ends. Arrays near an end are subtelomeric/terminal
    and are never chosen, but are reported. `ambiguous` = the runner-up internal array is at least ambiguity_ratio times as large as the chosen one."""
    arr = [dict(start=a, end=b, size=c, span=b - a, n_features=n, end_distance=min(a, length - b)) for a, b, c, n in arrays if c >= min_array]
    for a in arr:
        a['internal'] = a['start'] >= min_end_distance and (length - a['end']) >= min_end_distance
    internal = sorted([a for a in arr if a['internal']], key=lambda a: -a['size'])
    overall = max(arr, key=lambda a: a['size']) if arr else None
    best = internal[0] if internal else None
    second = internal[1]['size'] if len(internal) > 1 else 0
    return dict(best=best, best_overall=overall, n_arrays=len(arr), n_internal=len(internal), second_internal_size=second,
                ambiguous=bool(best and second >= ambiguity_ratio * best['size']), arrays=arr,
                overall_is_terminal=bool(best and overall and not overall['internal']))


def relation_to_inversion(c0, c1, i0, i1):
    """pericentric: centromere entirely inside the inversion; paracentric: entirely outside; otherwise a breakpoint lies in the centromere."""
    if any(pd.isna(x) for x in (c0, c1, i0, i1)):
        return ''
    if c1 < i0 or c0 > i1:
        return 'paracentric'
    if c0 >= i0 and c1 <= i1:
        return 'pericentric'
    return 'breakpoint within the centromere'
