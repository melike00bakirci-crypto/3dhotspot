"""Add regime + tier + caveat to the batch TSV, so every gene carries its own caveat."""
import csv

# --- B1 confound tanisali (DECISION-B1-DISORDER-CONFOUND-0001) -----------------
# b1_confound.py'nin null-vs-null olcumu: gercek kumelenme icermeyen, katli
# cekirdekten duzgun cekilmis bir kume ayni null'a karsi ne kadar yukselme
# uretiyor. ORAN TEK BASINA YANILTICI: gozlenen fazlasi minik olan genlerde
# patlar (PALB2 oran 41, ortusme 0.00). O yuzden ortusme ile kapili.
import json as _json
try:
    _B1 = {r["gene"]: r for r in _json.load(
        open("results_v3/b1_confound.json"))
        if r.get("explained") is not None}
except Exception:
    _B1 = {}

def _b1_status(g):
    r = _B1.get(g)
    if r is None:
        return "", "", "NOT_TESTED"
    e, o = r["explained"], r["overlap_median"]
    if e >= 0.50 and o >= 0.90:
        st = "CONFOUNDED"
    elif o >= 0.90:
        st = "BAND_CONFOUNDED"
    elif o < 0.50 and e < 0.30:
        st = "CLEAN"
    else:
        st = "PARTIAL"
    return round(e, 3), round(o, 3), st

IN = "results_v3/BATCH_DATA_ACMG_SF_v3.2_v3.tsv"
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None
def i(x):
    try: return int(float(x))
    except (TypeError, ValueError): return None

R = list(csv.DictReader(open(IN), delimiter="\t"))
FIELDS = list(R[0].keys()) + ["regime", "tier", "caveat",
                             "b1_explained", "b1_overlap", "b1_status"]

def classify(r):
    st, npv, nc = r["run_status"], i(r["N_P"]), i(r["n_centers"]) or 0
    rfp, cov = f(r["r_fp"]), f(r["fp_coverage"])
    fe, jc = f(r["fold_enr_at_r_hot"]), f(r["jaccard_median"])
    ol = str(r["binding_is_outlier"]).lower() == "true"

    if npv is None:
        return ("input_limited", "",
                f"Stage A could not assemble usable inputs: {r['stop_reason'][:90]}")
    if st == "UNDERPOWERED":
        return ("cannot_be_asked", "",
                f"N_P={npv}: the attainable p-floor exceeds q/m at any B (combinatorial "
                f"floor). The run says nothing about this gene.")
    if nc == 0:
        if st == "COMPLETED_NEGATIVE" and (r["stop_reason"] or "").startswith("INSUFF"):
            return ("input_limited", "",
                    "classified cohort is unanalysable (a binary class is empty)")
        return ("no_spatial_concentration", "",
                f"N_P={npv} with a passing power certificate and no significant centre: "
                f"the test could have rejected and did not.")
    if rfp is None:
        if (r["stop_reason"] or "") == "NO_ADMISSIBLE_FOOTPRINT":
            return ("refused_degenerate", "",
                    "r_contain = 0 — the centres ARE the pathogenic residues, so a "
                    "zero-radius union has no volume. Hotspot centres stand.")
        # [2026-10-30] Renamed from refused_excessive_coverage. "REFUSED" read as a
        # failure, which it is not: the centres stand and the finding is diffuse
        # pathogenicity. What is withheld is only the footprint GEOMETRY, and the
        # reason it carries no information is in the fold enrichment -- at fold ~ 1.0
        # the density inside the sphere equals the density outside, so a region
        # covering most of the protein adds nothing to "the variants are spread".
        fe_txt = f"fold enrichment {fe:.2f}" if fe is not None else "fold enrichment undefined"
        return ("diffuse_no_localised_region", "",
                f"pathogenicity is diffuse, not localised: the determined radius covers "
                f"{cov:.3f} of U_struct at {fe_txt}, so inside and outside the sphere have "
                f"the same pathogenic density. Hotspot centres STAND and are a finding; "
                f"only the footprint geometry is withheld, because a region covering most "
                f"of the protein is not a region.")

    score = sum([fe is not None and fe >= 1.5, jc is not None and jc >= 0.90,
                 not ol, cov is not None and cov < 0.30])
    tier = {4: "A", 3: "B", 2: "C", 1: "D", 0: "D"}[score]
    notes = []
    if fe is not None and fe < 1.05: notes.append("no enrichment over the labelled background")
    elif fe is not None and fe < 1.5: notes.append(f"weak enrichment ({fe:.2f})")
    if jc is not None and jc < 0.70: notes.append(f"fragile under centre loss (Jaccard {jc:.3f})")
    elif jc is not None and jc < 0.90: notes.append(f"moderate stability (Jaccard {jc:.3f})")
    if ol: notes.append("radius rests on a single outlying residue")
    if cov is not None and cov >= 0.30:
        notes.append(f"footprint covers {cov:.3f} of the structure, near the 0.50 limit")
    regime = "works_with_caveat" if tier == "D" else "works"
    return (regime, tier, "; ".join(notes) if notes else
            "meets all four quality criteria")

out = IN.replace(".tsv", "_classified.tsv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in R:
        reg, t, cav = classify(r)
        r["regime"], r["tier"], r["caveat"] = reg, t, cav
        r["b1_explained"], r["b1_overlap"], r["b1_status"] = _b1_status(r["gene"])
        w.writerow(r)
import collections

print(f"  {len(R)} gen -> {out}")
print("  rejim dagilimi:")
for k, v in collections.Counter(classify(r)[0] for r in R).most_common():
    print(f"    {k:<28} {v}")
