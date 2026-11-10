"""Build the analysis unit (one VCF = one species) and run the identity / design checks."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .findings import Findings, compress_ids, short_list


def norm_key(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _blank(x) -> bool:
    """True for None, NaN / pandas NA (an empty spreadsheet cell) and empty or whitespace-only strings."""
    if x is None:
        return True
    try:
        if pd.isna(x):
            return True
    except (TypeError, ValueError):
        pass
    return str(x).strip().lower() in ("", "nan")


def parse_sex(x) -> float:
    """Strict: only M or F (any case) are valid. Everything else -> NaN."""
    if _blank(x):
        return np.nan
    s = str(x).strip().upper()
    return 1.0 if s == "M" else 0.0 if s == "F" else np.nan


_UNKNOWN = {"", "?", "na", "nan", "none", "unknown", "unk", "-", "."}


def parse_parents(x, known=None) -> Tuple[str, Optional[str], Optional[str]]:
    """'sire,dam' -> ('ok', sire, dam); blank -> ('none',..).
    One parent unknown ('sire,', 'sire,?') -> ('partial', sire, None) / ('partial', None, dam).
    A single ID that is a known sample -> ('partial_single', id, None); its sex decides sire or dam.
    Anything else -> ('malformed',..)."""
    if _blank(x):
        return "none", None, None
    toks = [t.strip() for t in str(x).split(",")]
    if len(toks) == 2:
        a, b = [None if t.lower() in _UNKNOWN else t for t in toks]
        if a and b:
            return "ok", a, b
        if a or b:
            return "partial", a, b
        return "none", None, None
    if len(toks) == 1 and known is not None and toks[0] in known:
        return "partial_single", toks[0], None
    return "malformed", None, None


@dataclass
class Family:
    id: str
    sire: int                  # sample index (column in the VCF sample order)
    dam: int
    offspring: np.ndarray
    sex: np.ndarray            # 1 = male, 0 = female


@dataclass
class Unit:
    names: List[str]
    sex: np.ndarray
    species: str
    strata: np.ndarray                      # object labels: 'fam:<sire>x<dam>', 'pool', 'parent'
    contrast_idx: np.ndarray                # samples used by the sex contrasts (parents excluded)
    parent_idx: np.ndarray
    families: List[Family]                  # families whose transmission test can run
    component: Dict[str, int]               # stratum label -> independent-family component id
    rows: pd.DataFrame                      # sample-table rows aligned with names
    ecogroup: Optional[str] = None
    suggestions: List[Tuple[str, str, str]] = field(default_factory=list)   # (sample, listed parent, vcf name)
    excluded_idx: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))   # left out of the sex contrasts
    pairs: Dict[str, Tuple[Optional[str], Optional[str]]] = field(default_factory=dict)   # brood label -> (sire, dam)
    shared_parents: Dict[str, List[str]] = field(default_factory=dict)   # parent -> broods that share it


def _components(pairs: Dict[str, Tuple[Optional[str], Optional[str]]]) -> Dict[str, int]:
    """Union-find over parent names: families that share a (known) parent are not independent replicates."""
    parent = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for sire, dam in pairs.values():
        ps = [p for p in (sire, dam) if p]
        for p in ps:
            find(p)
        if len(ps) == 2:
            parent[find(ps[0])] = find(ps[1])
    roots = {}
    return {lab: roots.setdefault(find([p for p in (s, d) if p][0]), len(roots)) for lab, (s, d) in pairs.items()}


def pedigree_summary(names, rows, F: Findings) -> None:
    """Always reported (even when the run is blocked): how complete is the pedigree information?"""
    n_ok = n_none = n_mal = n_part = 0
    pairs, listed = {}, set()
    known = set(rows.index)
    for s in names:
        status, sire, dam = parse_parents(rows.at[s, "Parents"] if "Parents" in rows else None, known)
        if status in ("partial", "partial_single"):
            n_part += 1
            listed.update(p for p in (sire, dam) if p)
        elif status == "ok":
            n_ok += 1
            pairs[(sire, dam)] = pairs.get((sire, dam), 0) + 1
            listed.update([sire, dam])
        elif status == "malformed":
            n_mal += 1
        else:
            n_none += 1
    in_vcf = len(listed & set(names))
    if n_ok:
        msg = (f"Pedigree information in the sample table: {n_ok} of {len(names)} samples list 'sire,dam' "
               f"({len(pairs)} distinct parent pair{'s' if len(pairs) != 1 else ''}: "
               + "; ".join(f"{a} x {b} ({n} offspring)" for (a, b), n in list(pairs.items())[:6])
               + (" ..." if len(pairs) > 6 else "") + f"; {in_vcf} of the {len(listed)} listed parents are in this VCF); "
               f"{n_none} list no parents" + (f"; {n_part} list only one parent" if n_part else "")
               + (f"; {n_mal} are malformed" if n_mal else "") + ".")
    else:
        msg = (f"No sample in this VCF lists parents in the sample table ({n_none} blank"
               + (f", {n_mal} malformed" if n_mal else "") + "): all animals are treated as unrelated, so only the "
               "sex contrasts will run and the transmission (pedigree) test is off. Enter 'sire,dam' in the "
               "Parents column of the offspring (both parents must also be in this VCF) to enable it.")
    F.info("pedigree_summary", msg)


def _shared(pairs) -> Dict[str, List[str]]:
    users: Dict[str, List[str]] = {}
    for lab, (a, b) in pairs.items():
        for p_ in (a, b):
            if p_:
                users.setdefault(p_, []).append(lab)
    return {p_: v for p_, v in users.items() if len(v) > 1}


def build_unit(vcf_samples: List[str], tables, cfg, F: Findings) -> Optional[Unit]:
    st = tables.sample
    dup = st[st.SampleID.duplicated(keep=False) & st.SampleID.isin(vcf_samples)]
    if len(dup):
        ids = sorted(dup.SampleID.unique())
        F.block("duplicate_table_rows", f"{len(ids)} samples appear on more than one row of the sample table: "
                + compress_ids(ids), ids,
                why="With two rows it is unclear which sex/species/parents belong to the fish.",
                fix="Keep one row per SampleID in the SampleDatabase.",
                rows=[dict(sample=i, column="SampleID", current="duplicated", needed="one row") for i in ids])
    rows = st.drop_duplicates("SampleID").set_index("SampleID")

    missing = [s for s in vcf_samples if s not in rows.index]
    if missing:
        keys = {}
        for tid in rows.index:
            keys.setdefault(norm_key(tid), []).append(tid)
        near = {m: keys[norm_key(m)][0] for m in missing if len(keys.get(norm_key(m), [])) == 1}
        hint = (" Likely matches (same letters/digits, different case or punctuation): "
                + "; ".join(f"{m} -> {t}" for m, t in list(near.items())[:6]) + ("; ..." if len(near) > 6 else "") + ".") if near else ""
        F.block("sample_not_in_table", f"{len(missing)} of {len(vcf_samples)} samples in the VCF have no row in the "
                f"sample table: {compress_ids(missing)}.{hint}", missing,
                why="VCF sample names are matched to the SampleDatabase 'SampleID' column exactly; without a row "
                    "there is no sex, species or parentage for the fish.",
                fix="Add these samples to the SampleDatabase, or correct the ID so it matches the VCF header "
                    "(to rename on the VCF side use bcftools reheader).",
                rows=[dict(sample=m, column="SampleID", current="(no row)", needed=near.get(m, "a matching SampleID row"))
                      for m in missing])
    present = [s for s in vcf_samples if s in rows.index]
    r = rows.loc[present]
    if len(present):
        pedigree_summary(present, r, F)

    # ---- sex: every sample must be exactly M or F
    sex = r["Sex"].map(parse_sex)
    bad = list(r.index[sex.isna()])
    if bad:
        vals = r.loc[bad, "Sex"].map(lambda x: "blank" if _blank(x) else f"'{x}'").value_counts()
        breakdown = ", ".join(f"{k}: {v}" for k, v in vals.items())
        have_parents = sum(not _blank(r.at[s, "Parents"]) for s in bad) if "Parents" in r else 0
        F.block("sex_not_mf", f"{len(bad)} of {len(present)} samples ({len(bad) / len(present):.0%}) have no valid sex in "
                f"the sample table ({breakdown}). Affected: {compress_ids(bad)}.", bad,
                why="Sex is the trait being mapped. Samples whose sex is blank, U or uncertain (F?, M?) cannot be tested, and "
                    "the rule is that every sample in the VCF must be exactly M or F.",
                fix="In the SampleDatabase 'Sex' column enter M or F for each affected SampleID (the full list is in "
                    "fix_list.tsv). A fish whose sex cannot be determined has to be removed from the VCF."
                    + (f" ({have_parents} of them already list parents.)" if have_parents else ""),
                rows=[dict(sample=s, column="Sex", current="" if _blank(r.at[s, "Sex"]) else str(r.at[s, "Sex"]),
                           needed="M or F") for s in bad])

    # ---- species: all identical
    sp = r["Species"].map(lambda x: "" if _blank(x) else " ".join(str(x).split()))
    if (sp == "").any():
        b = list(sp.index[sp == ""])
        F.block("species_missing", f"{len(b)} of {len(present)} samples have no Species: {compress_ids(b)}.", b,
                why="The species check requires every sample to name its species.",
                fix="Fill the SampleDatabase 'Species' column for these samples (full list in fix_list.tsv).",
                rows=[dict(sample=s, column="Species", current="", needed="the species name") for s in b])
    distinct = sorted({s.lower() for s in sp if s})
    species_name = next((s for s in sp if s), "")
    if len(distinct) > 1:
        counts = sp[sp != ""].value_counts()
        desc = "; ".join(f"'{k}': {v} sample{'s' if v != 1 else ''} ({compress_ids(list(sp.index[sp == k]), 4)})"
                         for k, v in counts.items())
        F.block("species_mismatch", f"This VCF contains {len(counts)} different Species values: {desc}.", list(sp.index),
                why="Sex-determination systems and population structure differ between species, so one VCF must hold "
                    "a single species (parents included; cross families are not analysed).",
                fix="If these are spelling variants of one species, make the Species column identical for all of them. "
                    "Otherwise split the VCF by species (e.g. bcftools view -S <ids.txt>) and analyse each file.",
                rows=[dict(sample=s, column="Species", current=sp[s], needed="one species per VCF")
                      for s in sp.index if sp[s] != counts.index[0]])
    if F.blocked:
        return None

    names = list(vcf_samples)
    n = len(names)
    col = {s: i for i, s in enumerate(names)}
    key_full = {}
    key_nosuffix = {}
    for s in names:
        k = norm_key(s)
        key_full.setdefault(k, []).append(s)
        key_nosuffix.setdefault(re.sub(r"[mf]$", "", k), []).append(s)
    collide = [v for v in key_full.values() if len(v) > 1]
    if collide:
        F.warn("possible_duplicate_ids", "Sample names that differ only by case/punctuation (same fish "
               "sequenced twice?): " + "; ".join("/".join(v) for v in collide[:6]), [x for v in collide for x in v])

    sexv = r["Sex"].map(parse_sex).reindex(names).to_numpy(float)
    suggestions: List[Tuple[str, str, str]] = []
    fam_members: Dict[Tuple[Optional[str], Optional[str]], List[int]] = {}
    malformed, unresolved, partial = [], {}, {}
    known_ids = set(rows.index)
    for s in names:
        status, sire, dam = parse_parents(r.at[s, "Parents"] if "Parents" in r else None, known_ids)
        if status == "partial_single":                     # one ID given: its sex says whether it is the sire or the dam
            sx = parse_sex(rows.at[sire, "Sex"]) if sire in rows.index else np.nan
            if sx == 1.0:
                status = "partial"
            elif sx == 0.0:
                status, sire, dam = "partial", None, sire
            else:
                status = "malformed"
        if status == "malformed":
            malformed.append(f"{s} ({r.at[s, 'Parents']})")
        elif status in ("ok", "partial"):
            fam_members.setdefault((sire, dam), []).append(col[s])
    if malformed:
        F.warn("parents_malformed", "Parents field is not 'sire,dam' (treated as no pedigree): "
               + short_list(malformed), [m.split(" ")[0] for m in malformed])

    strata = np.array(["pool"] * n, dtype=object)
    families: List[Family] = []
    pairs: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    parent_set = set()
    for (sire, dam), idx in sorted(fam_members.items(), key=lambda kv: (kv[0][0] or "", kv[0][1] or "")):
        label = f"fam:{sire or '?'}x{dam or '?'}"
        pairs[label] = (sire, dam)
        for i in idx:
            strata[i] = label                               # full (or half) sibs stay together in the sex contrasts
        listed = [p for p in (sire, dam) if p]
        # a parent that IS in the VCF is a relative of the offspring, never an unrelated animal in the contrasts
        parent_set.update(col[p] for p in listed if p in col)
        if sire is None or dam is None:
            partial[label] = (listed[0], "dam" if dam is None else "sire", idx)
            continue
        missing_p = [p for p in (sire, dam) if p not in col]
        for p in missing_p:
            hit = key_full.get(norm_key(p)) or key_nosuffix.get(re.sub(r"[mf]$", "", norm_key(p)))
            if hit and len(hit) == 1:
                suggestions.append((label, p, hit[0]))
        if missing_p:
            unresolved[label] = (missing_p, len(idx), [("sire" if p == sire else "dam") for p in missing_p])
            continue
        # both parents are in the VCF: check their recorded sex
        s_sex, d_sex = sexv[col[sire]], sexv[col[dam]]
        if not (s_sex == 1.0 and d_sex == 0.0):
            F.warn("parent_sex_mismatch", f"{label}: first parent should be male (sire) and second female (dam); "
                   f"table says sire={'M' if s_sex == 1 else 'F'}, dam={'M' if d_sex == 1 else 'F'}. "
                   "Transmission test skipped for this family.", [sire, dam])
            continue
        off = np.array(idx)
        families.append(Family(label, col[sire], col[dam], off, sexv[off]))
    if unresolved:
        n_off = sum(v[1] for v in unresolved.values())
        desc = "; ".join(f"{k.replace('fam:', '')}: {v[1]} offspring, {'/'.join(v[2])} {', '.join(v[0])} missing"
                         for k, v in unresolved.items())
        F.warn("parent_not_in_vcf", f"{len(unresolved)} brood(s) list a parent that is not in this VCF ({n_off} offspring): "
               f"{desc}.", [p for v in unresolved.values() for p in v[0]],
               why="The transmission test needs both parents' genotypes to find informative markers, so these broods are left "
                   "out of it. Their offspring are still analysed as full sibs in the sex contrasts (association, Fst, "
                   "heterozygosity) and each brood still counts as an independent unit there.",
               fix="If the missing parent was sequenced, add them to the VCF to include the brood in the transmission test; "
                   "otherwise nothing is needed.",
               rows=[dict(sample=p, column="VCF sample", current="not in VCF", needed="sequenced and in the VCF")
                     for v in unresolved.values() for p in v[0]])
    if partial:
        n_off = sum(len(v[2]) for v in partial.values())
        desc = "; ".join(f"{k.replace('fam:', '')}: {len(v[2])} offspring, {v[1]} unknown" for k, v in partial.items())
        F.warn("parent_unknown", f"{len(partial)} brood(s) list only one parent ({n_off} offspring): {desc}.",
               [names[i] for v in partial.values() for i in v[2]],
               why="The transmission test needs both parents, so these broods are left out of it. Offspring sharing the known "
                   "parent are kept together as one family in the sex contrasts instead of being mixed with unrelated animals.",
               fix="Enter both parents as 'sire,dam' in the Parents column (and sequence the missing one) to include the brood "
                   "in the transmission test.",
               rows=[dict(sample=names[i], column="Parents", current=str(r.at[names[i], "Parents"]), needed="sire,dam")
                     for v in partial.values() for i in v[2]])
    if suggestions:
        F.info("parent_id_suggestions", "Listed parents with no exact VCF match but a likely one (NOT used; "
               "fix the Parents column to use them): " + short_list(f"{p} -> {v}" for _, p, v in suggestions))

    for i in parent_set:
        strata[i] = "parent"
    parent_idx = np.array(sorted(parent_set), dtype=int)
    contrast_idx = np.array([i for i in range(n) if i not in parent_set], dtype=int)
    eco = None
    if "Ecogroup" in r:
        e = [str(x) for x in r["Ecogroup"].dropna().unique() if str(x).strip()]
        eco = e[0] if len(e) == 1 else None
    return Unit(names, sexv, species_name, strata, contrast_idx, parent_idx, families,
                _components(pairs) if cfg.family_independence == "parent" else {lab: i for i, lab in enumerate(pairs)},
                r.reindex(names), eco, suggestions, pairs=pairs, shared_parents=_shared(pairs))


def design_checks(unit: Unit, cfg, F: Findings) -> None:
    """Sample-number checks. Counts use the samples that enter the sex contrasts (parents excluded)."""
    sx = unit.sex[unit.contrast_idx]
    nm, nf = int((sx == 1).sum()), int((sx == 0).sum())
    if min(nm, nf) < cfg.min_per_sex_block:
        npar = len(unit.parent_idx)
        F.block("too_few_per_sex", f"Only {nm} males and {nf} females are available for the sex contrasts (at least "
                f"{cfg.min_per_sex_block} of each are needed" + (f"; the {npar} parents are used for the transmission "
                "test only and are not counted" if npar else "") + ").",
                why="With fewer than 3 animals of a sex no between-sex comparison is possible.",
                fix="Add sexed animals of the under-represented sex to the VCF (or call more offspring).")
    elif min(nm, nf) < cfg.min_per_sex_warn:
        F.warn("low_per_sex", f"{nm} males / {nf} females: power is low and permutation p-values are coarse "
               f"(>= {cfg.min_per_sex_good} per sex recommended).")
    labs = [l for l in np.unique(unit.strata[unit.contrast_idx]) if l.startswith("fam:")]
    for lab in labs:
        m = unit.strata == lab
        k = int(m.sum())
        if k < cfg.family_min_offspring_warn:
            F.warn("small_family", f"{lab} has only {k} offspring (>= {cfg.family_min_offspring_warn} recommended).")
    comps = {unit.component[l] for l in labs}
    usable = [f for f in unit.families
              if min(int((f.sex == 1).sum()), int((f.sex == 0).sum())) >= cfg.min_offspring_per_sex]
    if labs:
        n_ind = len(comps)
        if unit.shared_parents:
            who = "; ".join(f"{p_} is a parent of {len(v)} broods ({', '.join(x.replace('fam:', '') for x in v)})"
                            for p_, v in unit.shared_parents.items())
            if cfg.family_independence == "parent":
                F.info("shared_parents", f"{len(labs)} broods but only {n_ind} independent famil{'y' if n_ind == 1 else 'ies'}: {who}. "
                       "Broods that share a parent are counted as one family (family_independence: parent). If they really have "
                       "different parents, correct the Parents column of the sample sheet.")
            else:
                F.info("shared_parents", f"{who}. These broods are treated as independent (family_independence: pair): they are "
                       "separate meioses with different partners, but the shared parent's haplotypes are common to both, so agreement "
                       "between them shows consistency across partners, not that the locus is shared across unrelated parents.")
        if n_ind < cfg.min_independent_families:
            msg = f"{n_ind} independent famil{'y' if n_ind == 1 else 'ies'} (>= {cfg.min_independent_families} recommended)"
            msg += (f"; {len(usable)} of them can be used for the transmission test." if len(usable) < len(labs) else ".")
            if n_ind < cfg.min_support_units:
                msg += " Replication across families is not possible, so the top confidence tier is unreachable from pedigrees alone."
            F.warn("few_independent_families", msg)
    if labs and not usable:
        F.warn("no_usable_pedigree", "No family has both parents genotyped and enough offspring of each sex: the "
               "transmission test is off; only sex contrasts will be run.")
    unit.families = usable
