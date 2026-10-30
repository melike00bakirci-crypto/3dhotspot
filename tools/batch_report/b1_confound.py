"""B1 — null-vs-null testi.

Gercek kumelenme ICERMEYEN bir kalinti kumesi (katli cekirdekten duzgun
dagilimla cekilmis) ayni pozisyonel null'a karsi g(r) sismesi uretiyor mu?
Uretiyorsa gozlenen yukselme kumelenmeden degil katli-vs-duzensiz
kontrastindan geliyor olabilir.

SALT OKUNUR. Hicbir dosya degistirilmez, hicbir sonuc yeniden yazilmaz.
"""
import csv, glob, json, os, sys
import numpy as np
sys.path.insert(0, "src")
from sphare3d.spatial.ripley import build_grid, cohort_statistics, intensity

B_DIAG = 600          # tanisal: yukselmenin YERINI olcuyoruz, p-degeri degil
N_REPEAT = 5          # ayni gende bes bagimsiz cekim
PLDDT_FOLDED = 70.0
SHELL = 1.0
GRID_STEP = 1.0
SEED = 20261030

def load(d):
    U = list(csv.DictReader(open(f"{d}/03_STRUCTURE_QC/positional_universe.tsv"), delimiter="\t"))
    C = list(csv.DictReader(open(f"{d}/03_STRUCTURE_QC/classified_cohort.tsv"), delimiter="\t"))
    ok = [r for r in U if r["ca_usable"].strip().upper() in ("TRUE", "1")]
    if not ok:
        ok = U
    coords = np.array([[float(r["x_ca"]), float(r["y_ca"]), float(r["z_ca"])] for r in ok])
    plddt = np.array([float(r["plddt"]) for r in ok])
    plp = [r for r in C if r["class"].upper().startswith("P")]
    return coords, plddt, len(plp)

def observed_E(d):
    """Gercek kosunun E kumesi ve en buyuk null-ustu fazlasi."""
    p = f"{d}/04_GLOBAL_CLUSTERING/pair_correlation_plp.tsv"
    if not os.path.exists(p):
        return None
    T = list(csv.DictReader(open(p), delimiter="\t"))
    E, exc = [], 0.0
    for r in T:
        try:
            rad, go, p95 = float(r["radius_A"]), float(r["g_obs"]), float(r["g_null_p95"])
        except (TypeError, ValueError):
            continue
        if go > 1.0 and go > p95:
            E.append(rad)
            exc = max(exc, go - p95)
    return {"E": E, "n_E": len(E), "max_excess": exc,
            "lo": min(E) if E else None, "hi": max(E) if E else None}

def main():
    targets = sys.argv[1:] or None
    out = []
    for d in sorted(glob.glob("results_v3/*/*/FULL_RESULTS")):
        gene = d.split("/")[1]
        if os.path.basename(os.path.dirname(d)).endswith(".partial"):
            continue
        if targets and gene not in targets:
            continue
        if not os.path.exists(f"{d}/03_STRUCTURE_QC/positional_universe.tsv"):
            continue
        obs = observed_E(d)
        if obs is None:
            continue
        coords, plddt, n_plp = load(d)
        if n_plp < 3 or len(coords) < 50:
            continue
        folded = np.flatnonzero(plddt >= PLDDT_FOLDED)
        if len(folded) < n_plp:
            out.append(dict(gene=gene, skipped="katli altkume N_P'den kucuk"))
            continue
        grid, d_max, _ = build_grid(coords, GRID_STEP)
        lam, vol, _ = intensity(coords)
        reps = []
        for k in range(N_REPEAT):
            rng = np.random.default_rng(SEED + 1000 * k + abs(hash(gene)) % 997)
            pick = rng.choice(folded, size=n_plp, replace=False)
            stats = cohort_statistics(
                "PLP", coords[pick], coords, grid, vol, SHELL, B_DIAG,
                rng, f"b1_confound::{gene}::{k}")
            go, p95 = np.asarray(stats.g_obs), np.asarray(stats.g_null_p95)
            mask = (go > 1.0) & (go > p95)
            rs = grid[mask]
            reps.append(dict(n_E=int(mask.sum()),
                             lo=float(rs.min()) if rs.size else None,
                             hi=float(rs.max()) if rs.size else None,
                             max_excess=float(np.max(go - p95)) if go.size else 0.0,
                             E=[float(x) for x in rs]))
        ne = [r["n_E"] for r in reps]
        me = [r["max_excess"] for r in reps]
        # gozlenen E ile confound E'sinin ortusmesi
        obs_set = set(obs["E"])
        ov = [len(obs_set & set(r["E"])) / len(obs_set) if obs_set else 0.0 for r in reps]
        out.append(dict(
            gene=gene, M=len(coords), N_P=n_plp,
            folded_frac=float((plddt >= PLDDT_FOLDED).mean()),
            obs_n_E=obs["n_E"], obs_lo=obs["lo"], obs_hi=obs["hi"],
            obs_max_excess=round(obs["max_excess"], 4),
            conf_n_E_median=float(np.median(ne)), conf_n_E_min=min(ne), conf_n_E_max=max(ne),
            conf_max_excess_median=round(float(np.median(me)), 4),
            conf_max_excess_max=round(float(np.max(me)), 4),
            overlap_median=round(float(np.median(ov)), 3),
            explained=round(float(np.median(me)) / obs["max_excess"], 3) if obs["max_excess"] > 0 else None,
        ))
        r = out[-1]
        print(f"  {gene:9s} N_P={r['N_P']:3d} katli={r['folded_frac']:.2f} | "
              f"GOZLENEN n_E={r['obs_n_E']:3d} [{r['obs_lo']},{r['obs_hi']}] fazla={r['obs_max_excess']:7.3f} | "
              f"CONFOUND n_E={r['conf_n_E_median']:5.1f} fazla={r['conf_max_excess_median']:7.3f} | "
              f"ortusme={r['overlap_median']:.2f} aciklanan={r['explained']}", flush=True)
    json.dump(out, open("results_v3/b1_confound.json", "w"), indent=1)
    print(f"\n  {len(out)} gen islendi")

main()
