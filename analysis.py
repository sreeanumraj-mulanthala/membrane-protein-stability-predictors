#!/usr/bin/env python3
"""Recompute the reported agreement statistics from the released tables.

Measured values are read from the MutTMP Stab-pred files (alpha.csv, beta.csv)
and, optionally, from the published BamA tables (Tables S8 and S9 of
George et al.); neither is redistributed here. Seed 20261106.

    python analysis.py --muttmp-dir PATH --out FOLDER [--bama-measured FILE]
"""

import argparse
import csv
import math
import os
import random
import sys
import tempfile
from collections import OrderedDict, defaultdict
from decimal import Decimal, InvalidOperation

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 20261106
NBOOT = 10000
N_MIN_PROTEIN = 30
FOLDS = ("Alpha", "Beta")
ALPHA = 0.05

# key, printed name, column in predictions_corpus.csv, kcal/mol units
PREDICTORS = [
    ("ddgun", "DDGun", "DDGun_ddG", True),
    ("esm2", "ESM-2", "ESM-2_score", False),
    ("rasp", "RaSP", "RaSP_ddG", True),
    ("thermompnn", "ThermoMPNN", "ThermoMPNN_ddG", True),
]
THREE = ("ddgun", "esm2", "thermompnn")
COHORTS = [("LEGACY", "In_LEGACY", 271),
           ("PRIMARY", "In_PRIMARY", 246),
           ("PRIMARY_SV", "In_PRIMARY_SV", 242)]

# key, printed name, column in predictions_bama.csv, kcal/mol units
BAMA_METHODS = [
    ("ddgun", "DDGun", "DDGun_ddG", True),
    ("esm2", "ESM-2", "ESM-2_score", False),
    ("rasp", "RaSP", "RaSP_ddG", True),
    ("thermompnn", "ThermoMPNN", "ThermoMPNN_ddG", True),
    ("mcsm", "mCSM-membrane", "mCSM-membrane_ddG", True),
    ("muttmb", "MutTMB Stab-pred", "MutTMB_ddG", True),
]
PE, PG = "PC:PE", "PC:PG"
N_BAMA = 88


class Stop(Exception):
    pass


# ------------------------------------------------------------- statistics --
def midrank(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(a, b):
    n = len(a)
    if n < 3:
        return float("nan")
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    if da == 0 or db == 0:
        return float("nan")
    return num / (da * db)


def spearman(a, b):
    if len(a) < 3:
        return float("nan")
    return pearson(midrank(a), midrank(b))


def mae(pred, exp):
    return sum(abs(x - y) for x, y in zip(pred, exp)) / len(pred)


def isnan(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def index_interval(values):
    """95% percentile interval, taking the sorted values at positions
    int(0.025 n) and int(0.975 n). Used for the per-predictor intervals."""
    if len(values) < 100:
        return (float("nan"), float("nan"))
    v = sorted(values)
    return (v[int(0.025 * len(v))], v[int(0.975 * len(v))])


def pct(sorted_vals, p):
    """Percentile of a sorted list with linear interpolation."""
    if not sorted_vals:
        return float("nan")
    k = (len(sorted_vals) - 1) * p
    lo = int(k)
    hi = min(lo + 1, len(sorted_vals) - 1)
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (k - lo)


def interp_interval(series, family=1):
    """Percentile interval with interpolation; family > 1 gives the
    Bonferroni-adjusted interval. Used for the contrasts."""
    s = sorted(series)
    a = ALPHA / 2.0 / float(family)
    return pct(s, a), pct(s, 1.0 - a)


# ---------------------------------------------------- corpus: bootstraps --
def clustered_bootstrap(recs, nboot, rng, kcal):
    """Resample proteins with replacement, keeping each protein's records
    together. recs carry exp, pred and protein_key."""
    by = defaultdict(list)
    for r in recs:
        by[r["protein_key"]].append(r)
    keys = sorted(by)
    if len(keys) < 2:
        return {}
    sp, pe, ma = [], [], []
    for _ in range(nboot):
        drawn = []
        for _ in range(len(keys)):
            drawn.extend(by[keys[rng.randrange(len(keys))]])
        if len(drawn) < 3:
            continue
        e = [d["exp"] for d in drawn]
        p = [d["pred"] for d in drawn]
        s, q = spearman(e, p), pearson(e, p)
        if not math.isnan(s):
            sp.append(s)
        if not math.isnan(q):
            pe.append(q)
        if kcal:
            ma.append(mae(p, e))
    out = {"spearman": index_interval(sp), "pearson": index_interval(pe)}
    if kcal:
        out["mae"] = index_interval(ma)
    return out


def record_bootstrap(recs, nboot, rng):
    """Resample records within one protein."""
    n = len(recs)
    if n < 10:
        return (float("nan"), float("nan"))
    sp = []
    for _ in range(nboot):
        drawn = [recs[rng.randrange(n)] for _ in range(n)]
        v = spearman([d["exp"] for d in drawn], [d["pred"] for d in drawn])
        if not math.isnan(v):
            sp.append(v)
    if len(sp) < 100:
        return (float("nan"), float("nan"))
    sp.sort()
    return (sp[int(0.025 * len(sp))], sp[int(0.975 * len(sp))])


def by_protein(rows):
    out = {}
    for r in rows:
        out.setdefault(r["protein_key"], []).append(r)
    return out


def draw(groups, names, rng):
    samp = []
    for n in [rng.choice(names) for _ in names]:
        samp.extend(groups[n])
    return samp


def rho_of(sample, key):
    return spearman([r["preds"][key] for r in sample], [r["exp"] for r in sample])


def boot_pairs(rows, keys, nboot):
    """Every predictor scored on the same resampled proteins."""
    groups = by_protein(rows)
    names = sorted(groups)
    rng = random.Random(SEED)
    draws = []
    for _ in range(nboot):
        samp = draw(groups, names, rng)
        got = {k: rho_of(samp, k) for k in keys}
        if any(isnan(v) for v in got.values()):
            continue
        draws.append(got)
    return draws


def boot_folds(rows, keys, nboot):
    """Proteins resampled within each fold, the same draw for every
    predictor."""
    strata = {}
    for f in FOLDS:
        g = by_protein([r for r in rows if r["fold"] == f])
        strata[f] = (g, sorted(g))
    rng = random.Random(SEED)
    draws = []
    for _ in range(nboot):
        samp = {f: draw(strata[f][0], strata[f][1], rng) for f in FOLDS}
        got, bad = {}, False
        for k in keys:
            got[k] = {}
            for f in FOLDS:
                v = rho_of(samp[f], k)
                if isnan(v):
                    bad = True
                got[k][f] = v
        if not bad:
            draws.append(got)
    return draws


# --------------------------------------------------------------- loading --
def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    return [{(k or "").strip(): v for k, v in r.items()} for r in rows]


def fnum(s):
    s = (s or "").strip()
    if s == "" or s.lower() in ("na", "nan", "none", "null"):
        return None
    v = float(s)
    if not math.isfinite(v):
        raise Stop("not a finite number: %r" % s)
    return v


def unique_rows(rows, col, path):
    out = {}
    for r in rows:
        k = (r.get(col) or "").strip()
        if k in out:
            raise Stop("%s appears more than once in %s" % (k, path))
        out[k] = r
    return out


def load_muttmp(folder):
    """Measured values and fold, keyed by Mutation_ID = structure + '_' +
    Mutation_PDB, as in the annotation table."""
    out = {}
    for name in ("beta.csv", "alpha.csv"):
        path = os.path.join(folder, name)
        if not os.path.exists(path):
            raise Stop("not found: %s. Download alpha.csv and beta.csv from "
                       "the MutTMP Stab-pred website into one folder." % path)
        for r in read_csv(path):
            sid = (r.get("PDB ID") or "").strip()
            mut = (r.get("Mutation_PDB") or "").strip()
            fold = (r.get("Strcuture_Class") or r.get("Structure_Class") or "").strip()
            mid = sid + "_" + mut
            out.setdefault(mid, []).append(
                {"exp": fnum(r.get("change_delta_g_h2o")), "fold": fold})
    return out


def load_corpus(annotation, predictions, muttmp_dir):
    ann = read_csv(annotation)
    legacy = OrderedDict()
    for r in ann:
        if r.get("In_LEGACY") == "True":
            if r["Mutation_ID"] in legacy:
                raise Stop("%s appears more than once in %s"
                           % (r["Mutation_ID"], annotation))
            legacy[r["Mutation_ID"]] = r
    if len(legacy) != 271:
        raise Stop("the annotation table gives %d LEGACY records, expected 271"
                   % len(legacy))
    meas = load_muttmp(muttmp_dir)
    preds = unique_rows(read_csv(predictions), "Mutation_ID", predictions)
    recs = {}
    for mid, a in legacy.items():
        m = meas.get(mid, [])
        if len(m) != 1:
            raise Stop("%s occurs %d times in the MutTMP files, expected once"
                       % (mid, len(m)))
        if m[0]["exp"] is None or m[0]["fold"] not in FOLDS:
            raise Stop("%s: no measured value or no fold in the MutTMP files"
                       % mid)
        p = preds.get(mid)
        if p is None:
            raise Stop("%s is not in %s" % (mid, predictions))
        recs[mid] = {"Mutation_ID": mid, "exp": m[0]["exp"],
                     "fold": m[0]["fold"], "protein_key": a["UniProt"],
                     "protein": a["Protein_Name"],
                     "preds": {k: fnum(p.get(col)) for k, _l, col, _u in PREDICTORS},
                     "cohorts": {c: a.get(col) == "True" for c, col, _n in COHORTS}}
    return recs


# ---------------------------------------------------------- corpus: main --
def corpus_cohort(recs, cohort, nboot):
    members = sorted(m for m, r in recs.items() if r["cohorts"][cohort])
    common = sorted(m for m in members
                    if all(recs[m]["preds"][k] is not None
                           for k, _l, _c, _u in PREDICTORS))

    def build(key, ids):
        return [{"exp": recs[m]["exp"], "pred": recs[m]["preds"][key],
                 "fold": recs[m]["fold"], "protein_key": recs[m]["protein_key"]}
                for m in sorted(ids)]

    rng = random.Random(SEED)
    metrics, lopo = [], []

    def row(label, level, name, rs, kcal, ci, resampling):
        e = [r["exp"] for r in rs]
        p = [r["pred"] for r in rs]
        out = OrderedDict(cohort=cohort, predictor=label, level=level,
                          name=name, n=len(rs), spearman=spearman(e, p),
                          spearman_lo=ci.get("spearman", ("", ""))[0],
                          spearman_hi=ci.get("spearman", ("", ""))[1],
                          pearson=pearson(e, p) if kcal else "",
                          pearson_lo=ci.get("pearson", ("", ""))[0] if kcal else "",
                          pearson_hi=ci.get("pearson", ("", ""))[1] if kcal else "",
                          mae=mae(p, e) if kcal else "",
                          mae_lo=ci.get("mae", ("", ""))[0] if kcal else "",
                          mae_hi=ci.get("mae", ("", ""))[1] if kcal else "",
                          resampling=resampling)
        metrics.append(out)

    head = {}
    for key, label, _c, kcal in PREDICTORS:
        rs = build(key, common)
        ci = clustered_bootstrap(rs, nboot, rng, kcal)
        row(label, "pooled", "all", rs, kcal, ci, "proteins")
        head[key] = spearman([r["exp"] for r in rs], [r["pred"] for r in rs])
    for key, label, _c, kcal in PREDICTORS:
        rs = build(key, common)
        for fold in sorted({r["fold"] for r in rs}):
            s = [r for r in rs if r["fold"] == fold]
            ci = clustered_bootstrap(s, nboot, rng, kcal)
            row(label, "fold", fold, s, kcal, ci, "proteins within fold")

    byprot = defaultdict(list)
    for m in common:
        r = recs[m]
        byprot[(r["protein_key"], r["protein"], r["fold"])].append(m)
    order = sorted(byprot.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for (up, name, fold), ids in order:
        if len(ids) < N_MIN_PROTEIN:
            continue
        for key, label, _c, kcal in PREDICTORS:
            rs = build(key, ids)
            lo, hi = record_bootstrap(rs, nboot, rng)
            row(label, "protein", name, rs, kcal,
                {"spearman": (lo, hi)}, "records within protein")
    for (up, name, fold), ids in order:
        rest = sorted(set(common) - set(ids))
        if len(rest) < 3:
            continue
        for key, label, _c, _k in PREDICTORS:
            rs = build(key, rest)
            rho = spearman([r["exp"] for r in rs], [r["pred"] for r in rs])
            lopo.append(OrderedDict(cohort=cohort, predictor=label,
                                    removed=name, uniprot=up,
                                    n_removed=len(ids), n_left=len(rs),
                                    spearman=rho, change=rho - head[key]))

    contrasts = []
    rows = [recs[m] for m in common]
    for family, keys in (("three", [k for k, *_ in PREDICTORS if k in THREE]),
                         ("four", [k for k, *_ in PREDICTORS])):
        labels = {k: l for k, l, *_ in PREDICTORS}
        pooled = {k: rho_of(rows, k) for k in keys}
        draws = boot_pairs(rows, keys, nboot)
        npair = len(keys) * (len(keys) - 1) // 2
        for i, ka in enumerate(keys):
            for kb in keys[i + 1:]:
                series = [d[ka] - d[kb] for d in draws]
                lo, hi = interp_interval(series)
                blo, bhi = interp_interval(series, npair)
                contrasts.append(OrderedDict(
                    cohort=cohort, family=family, contrast="pair",
                    a=labels[ka], b=labels[kb],
                    difference=pooled[ka] - pooled[kb], lo=lo, hi=hi,
                    bonferroni_lo=blo, bonferroni_hi=bhi,
                    family_size=npair, replicates=len(draws)))
        fdraws = boot_folds(rows, keys, nboot)
        for k in keys:
            pa = rho_of([r for r in rows if r["fold"] == "Alpha"], k)
            pb = rho_of([r for r in rows if r["fold"] == "Beta"], k)
            series = [d[k]["Beta"] - d[k]["Alpha"] for d in fdraws]
            lo, hi = interp_interval(series)
            blo, bhi = interp_interval(series, len(keys))
            contrasts.append(OrderedDict(
                cohort=cohort, family=family,
                contrast="beta-barrel minus alpha-helical", a=labels[k], b="",
                difference=pb - pa, lo=lo, hi=hi, bonferroni_lo=blo,
                bonferroni_hi=bhi, family_size=len(keys),
                replicates=len(fdraws)))
    return metrics, lopo, contrasts, len(members), len(common)


# ------------------------------------------------------------------ BamA --
def load_bama(predictions, measured):
    preds = {}
    for v, r in unique_rows(read_csv(predictions), "Variant", predictions).items():
        preds[v] = {k: fnum(r.get(col)) for k, _l, col, _u in BAMA_METHODS}
    meas = {PE: {}, PG: {}}
    for r in read_csv(measured):
        v = (r.get("variant") or r.get("mutation") or "").strip()
        c = (r.get("condition") or "").strip()
        x = (r.get("ddG_F") or r.get("ddG_F_kcal_mol") or "").strip()
        if c not in meas:
            raise Stop("condition must be PC:PE or PC:PG, found %r" % c)
        if v in meas[c]:
            raise Stop("%s appears twice for %s" % (v, c))
        try:
            ok = Decimal(x).is_finite()
        except InvalidOperation:
            ok = False
        if not ok:
            raise Stop("%s %s: ddG_F must be a number, found %r" % (v, c, x))
        meas[c][v] = x
    variants = sorted(meas[PE])
    if set(meas[PE]) != set(meas[PG]) or len(variants) != N_BAMA:
        raise Stop("the measured file must give the same %d variants in both "
                   "conditions" % N_BAMA)
    missing = [v for v in variants if v not in preds]
    if missing:
        raise Stop("variants without predictions: %s" % missing[:5])
    for v in variants:
        for k, label, _c, _u in BAMA_METHODS:
            if preds[v][k] is None:
                raise Stop("no %s prediction for %s" % (label, v))
    return preds, meas, variants


def bama(preds, meas, variants, nboot):
    def cell_rng(key, target):
        return random.Random("%d|%s|%s" % (SEED, key, target))

    def boot(pairs, rng):
        n = len(pairs)
        if n < 10:
            return (float("nan"), float("nan"))
        sp = []
        for _ in range(nboot):
            d = [pairs[rng.randrange(n)] for _ in range(n)]
            v = spearman([x[0] for x in d], [x[1] for x in d])
            if not math.isnan(v):
                sp.append(v)
        if len(sp) < 100:
            return (float("nan"), float("nan"))
        sp.sort()
        return (sp[int(0.025 * len(sp))], sp[int(0.975 * len(sp))])

    out = []
    targets = [(PE, PE, [float(meas[PE][v]) for v in variants]),
               (PG, PG, [float(meas[PG][v]) for v in variants]),
               ("PC:PE minus PC:PG", "DIFFERENCE",
                [float(Decimal(meas[PE][v]) - Decimal(meas[PG][v]))
                 for v in variants])]
    for name, tag, exp in targets:
        for key, label, _c, kcal in BAMA_METHODS:
            pr = [preds[v][key] for v in variants]
            lo, hi = boot(list(zip(exp, pr)), cell_rng(key, tag))
            out.append(OrderedDict(
                target=name, method=label, n=len(variants),
                spearman=spearman(exp, pr), lo=lo, hi=hi,
                pearson=pearson(exp, pr),
                mae=(mae(pr, exp) if kcal and tag != "DIFFERENCE" else "")))
            if not kcal:
                out[-1]["pearson"] = ""
    diff = targets[2][2]
    summary = OrderedDict(mean_difference=sum(diff) / len(diff),
                          positive=sum(1 for d in diff if d > 0),
                          n=len(diff))
    return out, summary


# ------------------------------------------------------------------- I/O --
def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(v) if isinstance(v, float) else v)
                        for k, v in r.items()})


def f3(x):
    return "n/a" if isnan(x) or x == "" else "%+.3f" % x


# --------------------------------------------------------------- self-test --
def selftest():
    results = []

    def check(name, ok):
        results.append(ok)
        print("%s  %s" % ("PASS" if ok else "FAIL", name))

    check("Spearman of a monotone pair is 1",
          abs(spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1) < 1e-12)
    check("Spearman of a reversed pair is -1",
          abs(spearman([1, 2, 3, 4], [4, 3, 2, 1]) + 1) < 1e-12)
    check("tied values share the mean rank",
          midrank([5.0, 1.0, 5.0, 2.0]) == [3.5, 1.0, 3.5, 2.0])
    check("a constant series gives an undefined correlation",
          math.isnan(spearman([1.0, 2.0, 3.0], [7.0, 7.0, 7.0])))
    check("interpolated percentile of 1..4 at 0.5 is 2.5",
          pct([1.0, 2.0, 3.0, 4.0], 0.5) == 2.5)
    ser = [float(i) for i in range(1000)]
    a, b = interp_interval(ser), interp_interval(ser, 3)
    check("the Bonferroni interval is wider", b[0] < a[0] and b[1] > a[1])
    check("a difference at printed precision has no rounding residue",
          float(Decimal("0.35") - Decimal("0.19")) == 0.16)
    recs = []
    for p, s in (("P1", 1.0), ("P2", -1.0), ("P3", 1.0)):
        for e in (-1.0, -2.0, -3.0, -4.0):
            recs.append({"exp": e, "pred": s * e, "fold": "Alpha",
                         "protein_key": p, "preds": {"a": s * e, "b": s * e}})
    c1 = clustered_bootstrap(recs, 300, random.Random(SEED), True)
    c2 = clustered_bootstrap(recs, 300, random.Random(SEED), True)
    check("the same seed gives the same interval", c1 == c2)
    d = boot_pairs(recs, ["a", "b"], 200)
    check("identical predictors give zero paired differences",
          all(x["a"] - x["b"] == 0 for x in d) and len(d) > 0)
    rows2 = []
    for p, f, sgn in (("A1", "Alpha", -1), ("A2", "Alpha", -1),
                      ("B1", "Beta", 1), ("B2", "Beta", 1)):
        for e in (-1.0, -2.0, -3.0, -4.0):
            rows2.append({"exp": e, "fold": f, "protein_key": p,
                          "preds": {"c": sgn * e}})
    fd = boot_folds(rows2, ["c"], 200)
    check("a known fold contrast of +2 is recovered",
          all(abs(x["c"]["Beta"] - x["c"]["Alpha"] - 2.0) < 1e-9 for x in fd))
    with tempfile.TemporaryDirectory() as tmp:
        def put(name, text):
            path = os.path.join(tmp, name)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            return path

        def stops(fn):
            try:
                fn()
            except Stop:
                return True
            return False

        cols = [c for _k, _l, c, _u in BAMA_METHODS]
        good = ",".join(["1.0"] * len(cols))
        vs = ["V%d" % i for i in range(N_BAMA)]
        pred_ok = put("p.csv", "Variant," + ",".join(cols) + "\n" +
                      "".join("%s,%s\n" % (v, good) for v in vs))
        meas_rows = "".join("%s,%s,0.5\n" % (v, c) for v in vs for c in (PE, PG))
        meas_ok = put("m.csv", "variant,condition,ddG_F\n" + meas_rows)
        check("clean BamA inputs are accepted",
              not stops(lambda: load_bama(pred_ok, meas_ok)))
        meas_nan = put("m2.csv", "variant,condition,ddG_F\n" +
                       meas_rows.replace("V0,PC:PE,0.5", "V0,PC:PE,NaN"))
        check("a NaN measured value stops the run",
              stops(lambda: load_bama(pred_ok, meas_nan)))
        pred_blank = put("p2.csv", open(pred_ok).read().replace(
            "V1," + good, "V1," + good[:-3]))
        check("a blank prediction stops the run",
              stops(lambda: load_bama(pred_blank, meas_ok)))
        pred_dup = put("p3.csv", open(pred_ok).read() + "V2," + good + "\n")
        check("a repeated identifier stops the run",
              stops(lambda: load_bama(pred_dup, meas_ok)))
    ok = all(results)
    print("SELFTEST %s %d of %d" % ("PASSED" if ok else "FAILED",
                                    sum(results), len(results)))
    return ok


def main():
    ap = argparse.ArgumentParser(description="Recompute the reported statistics")
    ap.add_argument("--muttmp-dir")
    ap.add_argument("--out")
    ap.add_argument("--annotation",
                    default=os.path.join(HERE, "data", "annotation_table.csv"))
    ap.add_argument("--predictions",
                    default=os.path.join(HERE, "data", "predictions_corpus.csv"))
    ap.add_argument("--bama-predictions",
                    default=os.path.join(HERE, "data", "predictions_bama.csv"))
    ap.add_argument("--bama-measured", default=None)
    ap.add_argument("--bootstrap", type=int, default=NBOOT)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if not a.muttmp_dir or not a.out:
        ap.error("--muttmp-dir and --out are required")
    if os.path.exists(a.out):
        sys.exit("STOP: %s exists; choose a new output folder." % a.out)
    try:
        recs = load_corpus(a.annotation, a.predictions, a.muttmp_dir)
        bama_in = (load_bama(a.bama_predictions, a.bama_measured)
                   if a.bama_measured else None)
    except Stop as e:
        sys.exit("STOP: %s" % e)

    os.makedirs(a.out)
    allm, alll, allc = [], [], []
    for cohort, _col, n_expected in COHORTS:
        m, l, c, n, nc = corpus_cohort(recs, cohort, a.bootstrap)
        if n != n_expected:
            sys.exit("STOP: %s has %d records, expected %d"
                     % (cohort, n, n_expected))
        allm += m
        alll += l
        allc += c
        print("%s: %d records, %d scored by all four predictors"
              % (cohort, n, nc))
        for r in m:
            if r["level"] == "pooled":
                print("  %-11s rho %s (%s to %s)"
                      % (r["predictor"], f3(r["spearman"]),
                         f3(r["spearman_lo"]), f3(r["spearman_hi"])))
    write_csv(os.path.join(a.out, "corpus_metrics.csv"), allm)
    write_csv(os.path.join(a.out, "corpus_lopo.csv"), alll)
    write_csv(os.path.join(a.out, "corpus_contrasts.csv"), allc)

    if bama_in:
        rows, summary = bama(*bama_in, nboot=a.bootstrap)
        write_csv(os.path.join(a.out, "bama_metrics.csv"), rows)
        print("BamA: %d variants; PC:PE minus PC:PG positive for %d, "
              "mean %+.4f kcal/mol" % (summary["n"], summary["positive"],
                                       summary["mean_difference"]))
        for r in rows:
            print("  %-18s %-16s rho %s (%s to %s)"
                  % (r["target"], r["method"], f3(r["spearman"]),
                     f3(r["lo"]), f3(r["hi"])))
    else:
        print("BamA: skipped (no --bama-measured file)")
    print("written to %s" % a.out)


if __name__ == "__main__":
    main()
