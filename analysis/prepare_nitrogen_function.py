"""Rebuild fixed functional scores offline from the released KO matrix."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/nitrogen_function"
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mass(members, matrix, samples):
    present = sorted(set(members) & set(matrix.index))
    return matrix.loc[present, samples].sum(axis=0) if present else pd.Series(0.0, index=samples)


def logged(values):
    values = np.asarray(values, float)
    positives = values[values > 0]
    delta = float(positives.min() / 2) if len(positives) else 1e-12
    return np.log(np.where(values > 0, values, delta)), delta, int((values == 0).sum())



def main():
    OUT.mkdir(parents=True, exist_ok=True)
    definition_path = ROOT / "data/function_definitions.json"
    cfg = json.loads(definition_path.read_text(encoding="utf-8"))
    meta = pd.read_csv(ROOT / "data_processed/metadata_confirmed.csv", dtype={"animal_id": str})
    samples = meta.sample_id.tolist()
    ko = pd.read_csv(ROOT / "data/KO_reconstructed_gene_relative.tsv", sep="\t", index_col=0)
    total_map = json.loads((ROOT / "data/gene_definition_evidence.json").read_text())["gene_relative"]["column_sums"]
    totals = pd.Series({s: total_map[s] for s in samples})
    ko_total = ko[samples].sum(axis=0)
    assert len(meta) == meta.animal_id.nunique() == 30
    assert meta.groupby("diet_group").size().to_dict() == dict.fromkeys("ABCDE", 6)
    assert np.isfinite(ko.to_numpy()).all() and (ko.to_numpy() >= 0).all()
    assert (ko_total > 0).all() and (ko_total <= totals * 1.000001).all()
    sets, carbon_sets = cfg["base_sets"], cfg["carbon_sets"]
    RATIOS = cfg["ratios"]
    for name in ["pathway_membership_summary.csv", "enzyme_membership.csv"]:
        pd.read_csv(ROOT / "data" / name).to_csv(OUT / name, index=False)
    data = meta[["sample_id", "animal_id", "diet_group", "initial_BW_kg"]].copy()
    data["all_gene_mass"] = totals.loc[samples].to_numpy()
    data["KO_gene_mass"] = ko_total.loc[samples].to_numpy()
    data["log_KO_annotation_load"] = np.log(ko_total / totals).to_numpy()
    desc, contributions, zero_rules, within = [], [], {}, data.iloc[:, :4].copy()
    for name, members in sets.items():
        v = mass(members, ko, samples)
        frac = v / totals
        data[name], delta, n_zero = logged(frac)
        within[name], within_delta, _ = logged(v / ko_total)
        data[name + "__fraction"] = frac.to_numpy()
        zero_rules[name] = {"all_gene_delta": delta, "KO_gene_delta": within_delta, "zero_n": n_zero}
        desc.append({
            "Feature_ID": name, "Catalog_KOs": len(members),
            "Present_KOs": len(set(members) & set(ko.index)),
            "Detected_animals": int((v > 0).sum()),
            "Mean_percent_all_genes": float(frac.mean() * 100),
            "Median_percent_all_genes": float(frac.median() * 100),
            "CV": float(frac.std() / frac.mean()) if frac.mean() > 0 else np.nan,
            "Zero_replacements": n_zero,
        })
        for k in members:
            vv = ko.loc[k, samples] if k in ko.index else pd.Series(0.0, index=samples)
            contributions.append({
                "Feature_ID": name, "KO": k, "Present_in_catalog": k in ko.index,
                "Detected_animals": int((vv > 0).sum()),
                "Mean_percent_all_genes": float((vv / totals).mean() * 100),
                "Cohort_set_mass_share": float(vv.sum() / v.sum()) if v.sum() > 0 else np.nan,
            })
    for name, (a, b) in RATIOS.items():
        data[name] = data[a] - data[b]
        within[name] = within[a] - within[b]
    features = list(sets) + list(RATIOS)
    assert len(features) == 23 and np.isfinite(data[features].to_numpy()).all()
    pairs, pair_scores = [], []
    for name, members in sets.items():
        for c_name, c_members in carbon_sets.items():
            shared = sorted(set(members) & set(c_members))
            n_members = sorted(set(members) - set(shared))
            c_exclusive = sorted(set(c_members) - set(shared))
            n_frac = mass(n_members, ko, samples) / totals
            c_frac = mass(c_exclusive, ko, samples) / totals
            yn, dn, nz = logged(n_frac)
            xc, dc, cz = logged(c_frac)
            pair_id = name + "__" + c_name
            pairs.append({
                "Pair_ID": pair_id, "Nitrogen_feature": name, "Carbon_feature": c_name,
                "Shared_KOs": "|".join(shared), "N_exclusive_KOs": "|".join(n_members),
                "C_exclusive_KOs": "|".join(c_exclusive),
                "N_positive_animals": int((n_frac > 0).sum()),
                "C_positive_animals": int((c_frac > 0).sum()),
                "N_zero_replacements": nz, "C_zero_replacements": cz,
                "N_zero_delta": dn, "C_zero_delta": dc,
                "Estimable": bool(np.std(yn) > 1e-12 and np.std(xc) > 1e-12),
            })
            for i, sample in enumerate(samples):
                pair_scores.append({
                    "Pair_ID": pair_id, "sample_id": sample, "animal_id": data.animal_id.iloc[i],
                    "diet_group": data.diet_group.iloc[i], "initial_BW_kg": data.initial_BW_kg.iloc[i],
                    "log_KO_annotation_load": data.log_KO_annotation_load.iloc[i],
                    "response": yn[i], "carbon": xc[i],
                })
    assert len(pairs) == 80 and len(pair_scores) == 2400
    for frame, name in [
        (data, "scores_per_animal.csv"), (within, "within_KO_scores_per_animal.csv"),
        (pd.DataFrame(desc), "functional_description.csv"),
        (pd.DataFrame(contributions), "KO_contribution.csv"),
        (pd.DataFrame(pairs), "carbon_pair_membership.csv"),
        (pd.DataFrame(pair_scores), "carbon_pair_scores.csv"),
    ]:
        frame.to_csv(OUT / name, index=False)

    cfg["prepared_files"] = {name: sha(OUT / name) for name in [
        "carbon_pair_membership.csv", "carbon_pair_scores.csv", "enzyme_membership.csv",
        "functional_description.csv", "KO_contribution.csv", "pathway_membership_summary.csv",
        "scores_per_animal.csv", "within_KO_scores_per_animal.csv",
    ]}
    cfg["inputs"] = {name: sha(ROOT / name) for name in [
        "data/function_definitions.json", "data/KO_reconstructed_gene_relative.tsv",
        "data/gene_definition_evidence.json", "data_processed/metadata_confirmed.csv",
        "data/enzyme_membership.csv", "data/pathway_membership_summary.csv",
    ]}
    cfg["historical_hashes"] = {}
    (OUT / "definition_freeze.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    print("PASS: 30 animals, 23 diet responses and 80 carbon-function pairs rebuilt offline.")

if __name__ == "__main__":
    main()
