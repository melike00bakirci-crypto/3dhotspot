"""Every number for every gene in the v3 batch, into one table."""
import csv, glob, json, os

def flat(p):
    try:
        j = json.load(open(p)); return j.get("payload") or j
    except Exception:
        return {}

wall = {}
for r in csv.DictReader(open("logs/acmg_v3_progress.tsv"), delimiter="\t"):
    wall[r["gene"]] = r.get("wall_s")

COLS = ["gene","run_status","outcome_type","wall_s","uniprot","mane_enst","mane_refseq",
        "M","N","N_P","N_B","n_conflict","D_max_A","r_hot","r_hot_runner_up","near_tie_margin",
        "n_centers","n_regions","coverage_at_r_hot","fold_enr_at_r_hot",
        "r_contain_exact","r_fp","fp_residues","fp_components","fp_volume_A3","fp_coverage",
        "binding_residue","binding_is_outlier","n_plp_outside_evidence",
        "jaccard_median","n_iterations","blocking","major","advisory","stop_reason"]

rows = []
# Gen basina TEK kosu: ayni gende birden fazla RUN_ID varsa (ornek: TSC1'in
# basarisiz koşusu + yeniden koşumu) en yeni zaman damgasini al. ".partial"
# dizinleri yarim kalmis kosulardir, hic alinmaz.
_by_gene = {}
for _d in glob.glob("results_v3/*/*/FULL_RESULTS"):
    _run = os.path.basename(os.path.dirname(_d))
    if _run.endswith(".partial"):
        continue
    _gene = _d.split("/")[1]
    if _gene not in _by_gene or _run > os.path.basename(os.path.dirname(_by_gene[_gene])):
        _by_gene[_gene] = _d

for d in [_by_gene[g] for g in sorted(_by_gene)]:
    g = d.split("/")[1]
    run = json.load(open(f"{d}/00_RUN_SUMMARY/run_summary.json"))
    coh = flat(f"{d}/02_CLINVAR/cohort_summary.json")
    rc = coh.get("residue_counts") or {}
    h1 = flat(f"{d}/03_STRUCTURE_QC/handoff_01.json")
    h2 = flat(f"{d}/06_FINAL_HOTSPOTS/handoff_02.json")
    h3 = flat(f"{d}/08_FINAL_FOOTPRINT/handoff_03.json")
    h4 = flat(f"{d}/09_ROBUSTNESS/handoff_04.json")
    try:
        dom = json.load(open(f"{d}/07_FOOTPRINT_RADIUS/footprint_domain.json"))
    except Exception:
        dom = {}
    cont = dom.get("containment") or {}
    prof = (h4.get("robustness_profile") or {})
    jac = (prof.get("jaccard_residues") or {}).get("median")
    fpcov = None
    try:
        dr = list(csv.DictReader(open(f"{d}/07_FOOTPRINT_RADIUS/footprint_descriptors.tsv"), delimiter="\t"))
        if dr: fpcov = float(dr[0]["coverage_fraction"])
    except Exception: pass
    sev = {"BLOCKING": 0, "MAJOR": 0, "ADVISORY": 0}
    for p in glob.glob(f"{d}/*/warnings_*.tsv"):
        try:
            for w in csv.DictReader(open(p), delimiter="\t"):
                k = (w.get("severity") or "").strip()
                if k in sev: sev[k] += 1
        except Exception: pass
    why = (run.get("negative_result") or {}).get("condition") or ""
    if not why and run.get("error"): why = str(run["error"]).split(":")[1][:70].strip() if ":" in str(run["error"]) else str(run["error"])[:70]
    if not why:
        for e in (run.get("escalations") or [])[:1]:
            why = str(e.get("ambiguity"))[:70]
    rows.append({
        "gene": g, "run_status": run.get("run_status"), "outcome_type": run.get("outcome_type"),
        "wall_s": wall.get(g, ""), "uniprot": coh.get("uniprot_acc") or h1.get("uniprot_acc"),
        "mane_enst": coh.get("mane_transcript"), "mane_refseq": coh.get("mane_refseq_nucleotide"),
        "M": h1.get("M"), "N": rc.get("N"), "N_P": rc.get("N_P"), "N_B": rc.get("N_B"),
        "n_conflict": rc.get("n_conflict"), "D_max_A": dom.get("D_max_A"),
        "r_hot": h2.get("hotspot_radius"), "r_hot_runner_up": h2.get("r_hot_runner_up_A"),
        "near_tie_margin": h2.get("r_hot_runner_up_margin_rel"),
        "n_centers": h2.get("n_significant_centers"), "n_regions": h2.get("n_hotspot_regions"),
        "coverage_at_r_hot": h2.get("coverage_fraction_at_r_hot"),
        "fold_enr_at_r_hot": h2.get("fold_enrichment_at_r_hot"),
        "r_contain_exact": cont.get("r_contain_exact_A"), "r_fp": h3.get("footprint_radius"),
        "fp_residues": h3.get("n_footprint_residues"), "fp_components": h3.get("n_components"),
        "fp_volume_A3": h3.get("volume_A3"), "fp_coverage": fpcov,
        "binding_residue": cont.get("binding_residue"),
        "binding_is_outlier": cont.get("binding_residue_is_outlier"),
        "n_plp_outside_evidence": cont.get("n_plp_outside_evidence"),
        "jaccard_median": jac, "n_iterations": (prof.get("iterations") or {}).get("completed"),
        "blocking": sev["BLOCKING"], "major": sev["MAJOR"], "advisory": sev["ADVISORY"],
        "stop_reason": why,
    })

out = "results_v3/BATCH_DATA_ACMG_SF_v3.2_v3.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in sorted(rows, key=lambda x: x["gene"]): w.writerow(r)
print(f"  {len(rows)} gen -> {out}")
print(f"  kolon: {len(COLS)}")
