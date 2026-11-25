#!/usr/bin/env python3
"""Refit the MutTMB Stab-pred model holding out one unit at a time (Figure 2).

Requires the MutTMB Stab-pred reproducibility package (--pkg), which is not
redistributed here. Paper values: Python 3.13.5, scikit-learn 1.8.0,
numpy 2.3.5, pandas 2.2.3.

    python heldout_ladder.py --pkg PATH --out ladder.csv
"""

import argparse
import json
import os
import platform
import sys

import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneGroupOut

EXPECTED = {"python": "3.13.5", "scikit-learn": "1.8.0", "numpy": "2.3.5",
            "pandas": "2.2.3"}

GROUPINGS = [("pair", "Pair_ID"), ("structure_file", "FILENAME"),
             ("pdb_entry", "_PDB_entry"), ("protein", "UniProt"),
             ("homology_cluster", "Homology_Group")]


def features_frame(df, features):
    return df.loc[:, features].apply(pd.to_numeric, errors="coerce")


def fit_predict(train, test, features, params, target):
    xtr = features_frame(train, features)
    xte = features_frame(test, features)
    med = xtr.median(numeric_only=True)
    xtr = xtr.fillna(med).fillna(0.0)
    xte = xte.fillna(med).fillna(0.0)
    ytr = pd.to_numeric(train[target], errors="raise").to_numpy(float)
    model = GradientBoostingRegressor(**params)
    model.fit(xtr, ytr)
    return model.predict(xte)


def held_out_by_group(df, groups, features, params, target):
    pred = np.full(len(df), np.nan)
    for tr, te in LeaveOneGroupOut().split(df, groups=groups):
        pred[te] = fit_predict(df.iloc[tr], df.iloc[te], features, params, target)
    return pred


def held_out_by_row(df, features, params, target):
    pred = np.full(len(df), np.nan)
    idx = np.arange(len(df))
    for i in idx:
        pred[i] = fit_predict(df.iloc[np.delete(idx, i)], df.iloc[[i]],
                              features, params, target)[0]
    return pred


def spearman(a, b):
    ra = pd.Series(np.asarray(a, float)).rank()
    rb = pd.Series(np.asarray(b, float)).rank()
    if len(ra) < 3 or ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def scores(pred, truth):
    return (float(np.corrcoef(pred, truth)[0, 1]),
            float(np.mean(np.abs(pred - truth))),
            spearman(pred, truth))


def environment():
    got = {"python": platform.python_version(),
           "scikit-learn": sklearn.__version__, "numpy": np.__version__,
           "pandas": pd.__version__}
    print("environment: " + ", ".join("%s %s" % kv for kv in got.items()))
    diff = [k for k in EXPECTED if got[k] != EXPECTED[k]]
    if diff:
        print("NOTE: differs from the environment used for the paper in %s "
              "(%s). Values may differ in the third decimal place."
              % (", ".join(diff),
                 ", ".join("%s %s" % (k, EXPECTED[k]) for k in diff)))


def main():
    ap = argparse.ArgumentParser(description="Held-out evaluation of the "
                                 "refitted MutTMB Stab-pred model")
    ap.add_argument("--pkg", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit("STOP: %s exists; choose a new output file." % a.out)
    for name in ("model_config.json", "selected_features.txt"):
        if not os.path.exists(os.path.join(a.pkg, name)):
            sys.exit("STOP: %s not found in %s" % (name, a.pkg))
    environment()

    with open(os.path.join(a.pkg, "model_config.json"), encoding="utf-8") as fh:
        cfg = json.load(fh)
    with open(os.path.join(a.pkg, "selected_features.txt"), encoding="utf-8") as fh:
        features = [x.strip() for x in fh.read().splitlines() if x.strip()]
    target = cfg.get("target_column", "change_delta_g_h2o")
    params = cfg["model_params"]
    df = pd.read_csv(os.path.join(a.pkg, "data", str(cfg["validation_training_file"])))
    df["_PDB_entry"] = df["FILENAME"].astype(str).str[:4]
    truth = pd.to_numeric(df[target], errors="coerce").to_numpy(float)
    direct = (df["Direction"].astype(str).str.strip().str.lower()
              == "direct").to_numpy()
    print("development table: %d rows, %d measured, %d features"
          % (len(df), int(direct.sum()), len(features)))

    preds = {}
    preds["in_sample"] = (fit_predict(df, df, features, params, target), 1)
    print("held out: row")
    preds["row"] = (held_out_by_row(df, features, params, target), len(df))
    for unit, col in GROUPINGS:
        print("held out: %s" % unit)
        g = df[col].astype(str).to_numpy()
        preds[unit] = (held_out_by_group(df, g, features, params, target),
                       int(len(pd.unique(g))))

    rows = []
    for unit, (p, ngroups) in preds.items():
        for basis, mask in (("all_rows", np.ones(len(df), bool)),
                            ("measured_only", direct)):
            r, mae, rho = scores(p[mask], truth[mask])
            rows.append({"held_out": unit, "scored": basis,
                         "rows": int(mask.sum()), "groups": ngroups,
                         "pearson_r": r, "mae": mae, "spearman": rho})
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    print()
    print(out.to_string(index=False, float_format=lambda x: "%.4f" % x))
    print("\nwritten to %s" % a.out)



if __name__ == "__main__":
    main()
