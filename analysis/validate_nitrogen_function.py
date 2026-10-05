"""Independently reconstruct scores and reproduce the fixed R statistics."""

from pathlib import Path
import csv
import hashlib
import json
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

from scipy import stats

OUT = ROOT / "results/nitrogen_function"
PRIOR = ROOT / "data"
cfg = json.loads((OUT / "definition_freeze.json").read_text(encoding="utf-8"))
errors = []


def check(actual, expected, label, tolerance=1e-7):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    assert np.array_equal(np.isnan(a), np.isnan(b)), label
    finite = np.isfinite(a) & np.isfinite(b)
    error = float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0
    errors.append(error)
    assert error < tolerance, (label, error)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    values = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    result = np.empty(len(p))
    result[order] = values
    return result


def fit(x, y):
    a = np.linalg.solve(x.T @ x, x.T)
    b = a @ y
    e = y - x @ b
    h = np.sum(x * a.T, axis=1)
    v = (a * (e / (1 - h))**2) @ a.T
    return b, v, e


def diet_fit(x, x0, y):
    b, v, e = fit(x, y)
    f = float(b[1:5] @ np.linalg.solve(v[1:5, 1:5], b[1:5]) / 4)
    r0 = y - x0 @ np.linalg.lstsq(x0, y, rcond=None)[0]
    r2 = (r0 @ r0 - e @ e)/(r0 @ r0)
    return b, v, f, r2


d = pd.read_csv(OUT / "scores_per_animal.csv", dtype={"animal_id": str})
w = pd.read_csv(OUT / "within_KO_scores_per_animal.csv")
meta = pd.read_csv(ROOT / "data_processed/metadata_confirmed.csv", dtype={"animal_id": str})
assert d.sample_id.tolist() == meta.sample_id.tolist()
assert d.animal_id.tolist() == meta.animal_id.tolist()
ko = pd.read_csv(PRIOR / "KO_reconstructed_gene_relative.tsv", sep="\t", index_col=0)
samples = d.sample_id.tolist()
gene_totals = json.loads((PRIOR / "gene_definition_evidence.json").read_text(
    encoding="utf-8"))["gene_relative"]["column_sums"]
total = np.array([gene_totals[s] for s in samples])
annotated = ko[samples].sum(axis=0).to_numpy()
check(d.all_gene_mass, total, "all-gene denominator", 1e-5)
check(d.KO_gene_mass, annotated, "KO-gene denominator", 1e-5)
check(d.log_KO_annotation_load, np.log(annotated/total), "annotation load")
# Verify the released, fixed EC/KO membership table.
enzyme_members = pd.read_csv(ROOT / "data/enzyme_membership.csv")
for name in cfg["EC_sets"]:
    expected = set(enzyme_members.loc[enzyme_members.Feature_ID == name, "KO"])
    assert expected == set(cfg["base_sets"][name])
mass_by_set = {}
for name, members in cfg["base_sets"].items():
    present = sorted(set(members) & set(ko.index))
    mass = ko.loc[present, samples].to_numpy().sum(axis=0)
    mass_by_set[name] = mass
    assert (mass > 0).all()
    check(d[name], np.log(mass/total), name + " score")
    check(w[name], np.log(mass/annotated), name + " within-KO score")
    check(d[name + "__fraction"], mass/total, name + " fraction")
for name, (a, b) in cfg["ratios"].items():
    check(d[name], np.log(mass_by_set[a]/mass_by_set[b]), name + " ratio")
    check(w[name], d[name], name + " denominator invariance")
features = cfg["features"]
bw = (d.initial_BW_kg.to_numpy() - d.initial_BW_kg.mean()) / 10
diet_dummies = np.column_stack([(d.diet_group == g).to_numpy(float) for g in "BCDE"])
x = np.column_stack([np.ones(30), diet_dummies, bw])
x0 = x[:, [0, 5]]
q0, q1 = np.linalg.qr(x0)[0], np.linalg.qr(x)[0]
orders = pd.read_csv(OUT / "diet_permutation_indices.csv").to_numpy(int) - 1
within_orders = pd.read_csv(OUT / "within_diet_permutation_indices.csv").to_numpy(int) - 1
assert orders.shape == within_orders.shape == (30, 9999)
assert np.all(np.sort(orders, axis=0) == np.arange(30)[:, None])
assert np.all(np.sort(within_orders, axis=0) == np.arange(30)[:, None])
assert np.all(d.diet_group.to_numpy()[within_orders] == d.diet_group.to_numpy()[:, None])

fits = {}
for filename, data in [("diet_omnibus.csv", d), ("diet_within_KO_sensitivity.csv", w)]:
    results = pd.read_csv(OUT / filename)
    assert results.Feature_ID.tolist() == features
    assert len(results) == 23
    for row in results.itertuples():
        assert row.Estimable
        y = data[row.Feature_ID].to_numpy()
        b, v, f, r2 = diet_fit(x, x0, y)
        check(f, row.F_HC3, row.Feature_ID + " HC3 F")
        check(stats.f.sf(f, 4, 24), row.P_HC3, row.Feature_ID + " HC3 P")
        check(r2, row.partial_R2, row.Feature_ID + " partial R2")
        nullfit = q0 @ (q0.T @ y)
        e0 = y - nullfit
        ystar = e0[orders] + nullfit[:, None]
        s0 = np.sum((ystar - q0 @ (q0.T @ ystar))**2, axis=0)
        s1 = np.sum((ystar - q1 @ (q1.T @ ystar))**2, axis=0)
        rss0, rss1 = e0 @ e0, np.sum((y-x@b)**2)
        observed = ((rss0-rss1)/4)/(rss1/24)
        permstat = ((s0-s1)/4)/(s1/24)
        p = (1 + np.sum(permstat >= observed - 1e-12))/10000
        check(p, row.P_permutation, row.Feature_ID + " residual permutation")
        if filename == "diet_omnibus.csv":
            fits[row.Feature_ID] = (b, v)
    check(bh(results.P_HC3), results.q_HC3, filename + " HC3 BH")
    check(bh(results.P_permutation), results.q_permutation, filename + " permutation BH")

contrasts = pd.read_csv(OUT / "diet_all_contrasts.csv")
assert len(contrasts) == 230
for row in contrasts.itertuples():
    reference = np.array([1] + [float(row.reference == g) for g in "BCDE"] + [0])
    comparison = np.array([1] + [float(row.comparison == g) for g in "BCDE"] + [0])
    c = comparison - reference
    b, v = fits[row.Feature_ID]
    beta, se = c @ b, np.sqrt(c @ v @ c)
    check([beta, se, np.exp(beta)], [row.log_ratio, row.SE_HC3, row.ratio], "diet contrast")
    check(2*stats.t.sf(abs(beta/se), 24), row.P, "contrast P")
    check(np.exp(beta + np.array([-1, 1])*stats.t.ppf(.975, 24)*se),
          [row.ratio_CI_low, row.ratio_CI_high], "contrast CI")
check(bh(contrasts.P), contrasts.q_all230, "230 contrast BH")
diet_loo = pd.read_csv(OUT / "diet_LOO.csv", dtype={"excluded_animal": str})
assert len(diet_loo) == 690
for row in diet_loo.itertuples():
    i = d.animal_id.tolist().index(row.excluded_animal)
    ids = np.arange(30) != i
    b, v, f, _ = diet_fit(x[ids], x0[ids], d[row.Feature_ID].to_numpy()[ids])
    check(f, row.F_HC3, "diet LOO F")
    check(stats.f.sf(f, 4, 23), row.P_HC3, "diet LOO P")
    check(b[1:5], [row.B_log_ratio, row.C_log_ratio, row.D_log_ratio, row.E_log_ratio],
          "diet LOO coefficients")

members = pd.read_csv(OUT / "carbon_pair_membership.csv", keep_default_na=False)
pair_data = pd.read_csv(OUT / "carbon_pair_scores.csv", dtype={"animal_id": str})
carbon = pd.read_csv(OUT / "carbon_function_associations.csv")
carbon_loo = pd.read_csv(OUT / "carbon_LOO.csv", dtype={"excluded_animal": str})
boot_ids = pd.read_csv(OUT / "bootstrap_indices.csv").to_numpy(int) - 1
assert len(carbon) == 80 and len(carbon_loo) == 2400
assert boot_ids.shape == (2000, 30)
assert all(np.all((d.diet_group.to_numpy()[boot_ids] == g).sum(axis=1) == 6) for g in "ABCDE")
for row in carbon.itertuples():
    m = members[members.Pair_ID == row.Pair_ID].iloc[0]
    n_set = set(cfg["base_sets"][m.Nitrogen_feature])
    c_set = set(cfg["carbon_sets"][m.Carbon_feature])
    overlap = n_set & c_set
    assert set(filter(None, m.Shared_KOs.split("|"))) == overlap
    n_only, c_only = sorted((n_set-overlap) & set(ko.index)), sorted((c_set-overlap) & set(ko.index))
    y = np.log(ko.loc[n_only, samples].to_numpy().sum(axis=0)/total)
    carbon_score = np.log(ko.loc[c_only, samples].to_numpy().sum(axis=0)/total)
    z = pair_data[pair_data.Pair_ID == row.Pair_ID]
    check(z.response, y, "exclusive nitrogen score")
    check(z.carbon, carbon_score, "exclusive carbon score")
    xx = np.column_stack([np.ones(30), carbon_score, diet_dummies, bw, np.log(annotated/total)])
    null = xx[:, [0, 2, 3, 4, 5, 6, 7]]
    xx_no_load = xx[:, :7]
    assert np.linalg.matrix_rank(xx) == 8 and row.Estimable
    b, v, e = fit(xx, y)
    beta, se = b[1], np.sqrt(v[1, 1])
    p = 2*stats.t.sf(abs(beta/se), 22)
    check([beta, se, p], [row.Beta, row.SE_HC3, row.P_HC3], "carbon HC3")
    check(beta + np.array([-1, 1])*stats.t.ppf(.975, 22)*se,
          [row.CI_lower, row.CI_upper], "carbon CI")
    check(beta*np.std(carbon_score, ddof=1)/np.std(y, ddof=1),
          row.Beta_standardized, "standardized carbon coefficient")
    qn, qx = np.linalg.qr(null)[0], np.linalg.qr(xx)[0]
    yr, xr = y-qn@(qn.T@y), carbon_score-qn@(qn.T@carbon_score)
    check(np.corrcoef(xr, yr)[0, 1], row.Partial_r, "residual partial r")
    bb, vv, _ = fit(xx_no_load, y)
    check([bb[1], 2*stats.t.sf(abs(bb[1]/np.sqrt(vv[1, 1])), 23)],
          [row.Beta_no_load, row.P_no_load], "no-load sensitivity")
    inv = np.linalg.inv(xx.T @ xx)
    aa = inv @ xx.T
    fitted0 = y - yr
    ystar = yr[within_orders] + fitted0[:, None]
    perm_b = aa[1] @ ystar
    rss = np.sum((ystar-qx@(qx.T@ystar))**2, axis=0)
    perm_t = perm_b/np.sqrt(rss/22*inv[1, 1])
    obs_t = beta/np.sqrt((e@e)/22*inv[1, 1])
    pperm = (1+np.sum(np.abs(perm_t) >= abs(obs_t)-1e-12))/10000
    check(pperm, row.P_permutation, "within-diet carbon permutation")
    loo = np.array([np.linalg.lstsq(np.delete(xx, i, axis=0), np.delete(y, i), rcond=None)[0][1]
                    for i in range(30)])
    saved_loo = carbon_loo[carbon_loo.Pair_ID == row.Pair_ID]
    assert saved_loo.excluded_animal.tolist() == d.animal_id.tolist()
    check(loo, saved_loo.Beta, "carbon LOO coefficients")
    check([np.mean(np.sign(loo) == np.sign(beta)), loo.min(), loo.max()],
          [row.LOO_sign_agreement, row.LOO_beta_min, row.LOO_beta_max], "carbon LOO summary")
    bootstrap = []
    for ids in boot_ids:
        bb, _, rank, _ = np.linalg.lstsq(xx[ids], y[ids], rcond=None)
        bootstrap.append(bb[1] if rank == 8 else np.nan)
    bootstrap = np.array(bootstrap)
    assert np.isfinite(bootstrap).sum() == row.Bootstrap_success
    check([np.mean(bootstrap[np.isfinite(bootstrap)] > 0),
           *np.nanquantile(bootstrap, [.025, .975])],
          [row.Bootstrap_positive_fraction, row.Bootstrap_CI_lower, row.Bootstrap_CI_upper],
          "carbon bootstrap", 1e-6)
for p, q in [("P_HC3", "q_HC3"), ("P_permutation", "q_permutation"), ("P_no_load", "q_no_load")]:
    check(bh(carbon[p]), carbon[q], "80-pair " + q)
assert np.array_equal(
    carbon.Supported_HC3_and_permutation.to_numpy(),
    ((carbon.q_HC3 < .05) & (carbon.q_permutation < .05)).to_numpy())
for mapping, prefix in [(cfg["inputs"], ROOT), (cfg["historical_hashes"], ROOT),
                        (cfg["prepared_files"], OUT)]:
    for name, expected in mapping.items():
        assert hashlib.sha256((prefix/name).read_bytes()).hexdigest() == expected, name
summary = {
    "status": "PASS", "diet_responses": 23, "diet_contrasts": 230,
    "diet_LOO_fits": 690, "carbon_pairs": 80, "carbon_LOO_fits": 2400,
    "carbon_bootstrap_fits": 160000, "permutations": 9999,
    "sources_and_scores_reconstructed": True, "all_HC3_and_P_reproduced": True,
    "all_BH_families_reproduced": True, "all_residual_permutations_reproduced": True,
    "all_LOO_and_bootstrap_summaries_reproduced": True,
    "historical_hashes_unchanged": True, "max_absolute_numeric_error": max(errors),
}
(OUT / "independent_validation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2), flush=True)
