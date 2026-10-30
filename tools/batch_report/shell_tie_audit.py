"""Salt okunur: kabuk seyrekligi + E kuralinin > vs >= duyarliligi."""
import csv, glob, json, math, os, statistics as st

ROOT = "results_v3"
R_FLOOR, R_CEIL, MARGIN = 5.0, 25.0, 1.0
TRIVIAL_R = 3.0
EPS = 1e-9
MISSING = "VERI_YOK"

def fl(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None

def jload(p):
    try:
        with open(p) as fh: return json.load(fh)
    except Exception:
        return None

def blocks(rs, step=1.0):
    """Ardisik farki <= step olan bitisik bloklar."""
    out = []
    for r in sorted(rs):
        if out and (r - out[-1][-1]) <= step + EPS:
            out[-1].append(r)
        else:
            out.append([r])
    return out

def longest(bl):
    """En uzun blok; esitlikte en kucuk baslangic (pipeline ile ayni kural)."""
    if not bl: return None
    best = bl[0]
    for b in bl[1:]:
        if (b[-1] - b[0]) > (best[-1] - best[0]):
            best = b
    return best

EMPTY_E   = "E_BOS"          # kume bos -> boru hattinda FT-1
EMPTY_CLIP = "KIRPMA_BOS"     # en uzun blok [5,25] disinda kaliyor

def domain_from(rs):
    """En uzun blok +/- MARGIN, [R_FLOOR, R_CEIL]'e kirpilmis.

    Donus: (lo,hi) | EMPTY_E (kume bos) | EMPTY_CLIP (kirpma sonrasi bos).
    """
    b = longest(blocks(rs))
    if b is None: return EMPTY_E
    lo = max(R_FLOOR, b[0] - MARGIN)
    hi = min(R_CEIL, b[-1] + MARGIN)
    if hi < lo: return EMPTY_CLIP
    return (round(lo, 6), round(hi, 6))

def fmt_dom(d):
    return d if isinstance(d, str) else f"[{d[0]:g},{d[1]:g}]"

def fmt_radii(rs):
    return ";".join(f"{r:g}" for r in rs) if rs else ""

def spearman(xs, ys):
    n = len(xs)
    if n < 3: return None
    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        rk = [0.0] * n
        i = 0
        while i < n:                      # bagli degerlere ortalama rank
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]: j += 1
            avg = (i + j) / 2.0
            for k in range(i, j + 1): rk[order[k]] = avg
            i = j + 1
        return rk
    rx, ry = rank(xs), rank(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return None if den == 0 else num / den

# ---- gen dizinlerini bul -------------------------------------------------
genes = {}
for d in sorted(glob.glob(f"{ROOT}/*/")):
    g = d.rstrip("/").split("/")[-1]
    # Yarim kalmis kosulari hic alma; ayni gende birden fazla tamamlanmis kosu
    # varsa (ornek: TSC1) RUN_ID zaman damgasiyla siralanir, en yenisi alinir.
    runs = sorted(r for r in glob.glob(f"{d}*/FULL_RESULTS")
                  if not os.path.basename(os.path.dirname(r)).endswith(".partial"))
    genes[g] = runs[-1] if runs else None

rows, notes = [], []
for g, fr in genes.items():
    rec = {"gene": g, "branch_taken": MISSING, "N_P": MISSING, "n_pairs": MISSING,
           "pcf_shell_width_A": MISSING, "n_shells": MISSING, "n_zero": MISSING,
           "n_tie": MISSING, "n_tie_trivial": MISSING, "n_tie_zero": MISSING,
           "n_tie_g_gt_1": MISSING, "tie_radii_g_gt_1": "",
           "zero_radii": "", "tie_radii": "",
           "tie_radii_trivial": "", "domain_gt": MISSING, "domain_ge": MISSING,
           "domain_changed": MISSING, "pivot_radii": "", "selected_r_hot": MISSING,
           "r_hot_in_domain_ge": MISSING, "missing_files": ""}
    if fr is None:
        rec["missing_files"] = "FULL_RESULTS"
        rows.append(rec); continue

    miss = []
    pcf = f"{fr}/04_GLOBAL_CLUSTERING/pair_correlation_plp.tsv"
    crd = f"{fr}/04_GLOBAL_CLUSTERING/candidate_radius_domain.json"
    rks = f"{fr}/04_GLOBAL_CLUSTERING/ripleys_k_summary.json"
    coh = f"{fr}/02_CLINVAR/cohort_summary.json"
    h02 = f"{fr}/06_FINAL_HOTSPOTS/handoff_02.json"

    # --- N_P / n_pairs
    c = jload(coh)
    if c is None: miss.append("cohort_summary.json")
    else:
        np_ = c.get("N_P")
        if np_ is not None:
            rec["N_P"] = int(np_)
            rec["n_pairs"] = int(np_) * (int(np_) - 1) // 2

    # --- kabuk genisligi
    k = jload(rks)
    if k is None: miss.append("ripleys_k_summary.json")
    else:
        w = fl(k.get("pcf_shell_width_A"))
        if w is not None: rec["pcf_shell_width_A"] = w

    # --- branch
    cd = jload(crd)
    if cd is None: miss.append("candidate_radius_domain.json")
    else: rec["branch_taken"] = cd.get("branch_taken") or MISSING

    # --- secilen r_hot
    h = jload(h02)
    if h is None: miss.append("handoff_02.json")
    else:
        rh = fl(h.get("hotspot_radius"))
        if rh is not None: rec["selected_r_hot"] = rh

    # --- A + B: pcf tablosu
    if not os.path.exists(pcf):
        miss.append("pair_correlation_plp.tsv")
    else:
        with open(pcf) as fh:
            T = list(csv.DictReader(fh, delimiter="\t"))
        rec["n_shells"] = len(T)
        zero, tie, tie_triv, e_gt, e_ge = [], [], [], [], []
        tie_zz, tie_gt1 = [], []
        for r in T:
            rad, go, p95 = fl(r["radius_A"]), fl(r["g_obs"]), fl(r["g_null_p95"])
            if rad is None or go is None or p95 is None: continue
            if go == 0.0: zero.append(rad)
            if abs(go - p95) <= EPS:
                if rad <= TRIVIAL_R + EPS:
                    tie_triv.append(rad)
                else:
                    tie.append(rad)
                    (tie_gt1 if go > 1.0 else tie_zz).append(rad)
            if go > 1.0:
                if go > p95: e_gt.append(rad)
                if go >= p95 - EPS: e_ge.append(rad)
        rec["n_zero"], rec["n_tie"], rec["n_tie_trivial"] = len(zero), len(tie), len(tie_triv)
        rec["n_tie_zero"] = len(tie_zz)
        rec["n_tie_g_gt_1"] = len(tie_gt1)
        rec["tie_radii_g_gt_1"] = fmt_radii(tie_gt1)
        rec["zero_radii"], rec["tie_radii"] = fmt_radii(zero), fmt_radii(tie)
        rec["tie_radii_trivial"] = fmt_radii(tie_triv)
        dg, de = domain_from(e_gt), domain_from(e_ge)
        rec["domain_gt"], rec["domain_ge"] = fmt_dom(dg), fmt_dom(de)
        rec["domain_changed"] = "TRUE" if dg != de else "FALSE"
        # pivot: E_ge'ye giren ama E_gt'de olmayan yaricaplar
        rec["pivot_radii"] = fmt_radii(sorted(set(e_ge) - set(e_gt)))
        rh = rec["selected_r_hot"]
        if not isinstance(rh, float):
            rec["r_hot_in_domain_ge"] = MISSING
        elif isinstance(de, str):
            rec["r_hot_in_domain_ge"] = f"TANIMSIZ({de})"
        else:
            rec["r_hot_in_domain_ge"] = "TRUE" if de[0] - EPS <= rh <= de[1] + EPS else "FALSE"

    rec["missing_files"] = ";".join(miss)
    rows.append(rec)

COLS = ["gene", "branch_taken", "N_P", "n_pairs", "pcf_shell_width_A", "n_shells",
        "n_zero", "n_tie", "n_tie_trivial", "n_tie_zero", "n_tie_g_gt_1",
        "zero_radii", "tie_radii", "tie_radii_trivial", "tie_radii_g_gt_1",
        "domain_gt", "domain_ge", "domain_changed",
        "pivot_radii", "selected_r_hot", "r_hot_in_domain_ge", "missing_files"]
out_tsv = f"{ROOT}/SHELL_TIE_AUDIT_v3.tsv"
with open(out_tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows: w.writerow(r)

# ---- OZET ---------------------------------------------------------------
L = []
A = L.append
have = [r for r in rows if r["n_shells"] != MISSING]
nodata = [r for r in rows if r["n_shells"] == MISSING]
pcfd = [r for r in have if r["branch_taken"] == "pcf_derived"]
fb = [r for r in have if r["branch_taken"] == "fallback_envelope"]
tie1 = [r for r in have if isinstance(r["n_tie"], int) and r["n_tie"] >= 1]
ch = [r for r in have if r["domain_changed"] == "TRUE"]
ch_pcf = [r for r in ch if r["branch_taken"] == "pcf_derived"]

A("# Kabuk seyrekligi ve E kurali esitlik duyarliligi")
A("")
A(f"Kok: `{ROOT}/` · {len(rows)} gen dizini · PCF tablosu okunan {len(have)} · "
  f"{MISSING} {len(nodata)}")
A("")
A("E kriteri karsilastirmasi: `E_gt = g_obs > 1 AND g_obs > g_null_p95` (mevcut kural) "
  "karsisinda `E_ge = g_obs > 1 AND g_obs >= g_null_p95`. Her iki kumede bitisik bloklar "
  "(ardisik fark <= 1 A), en uzun blok (esitlikte en kucuk baslangic), +/-1 A marj, "
  "[5,25] kirpma. Not: boru hattinin kendi ust kirpmasi `min(25, 0.25*D_max)`; burada "
  "istenen sabit [5,25] kullanildi.")
A("")
A("## Esitlik sayilari")
A("")
A("`g_obs > 1` sarti nedeniyle **yalnizca `g_obs > 1` olan esitlikler E kumesini "
  "etkileyebilir**; `g_obs == 0` esitlikleri (bos kabuk, null p95 de 0) E'ye hic giremez.")
A("")
tie_eff = [r for r in have if isinstance(r["n_tie_g_gt_1"], int) and r["n_tie_g_gt_1"] >= 1]
A(f"- `n_tie >= 1` (r > 3 A, trivial haric): **{len(tie1)}** / {len(have)} gen")
A(f"- bunlardan `g_obs > 1` olan esitligi bulunan: **{len(tie_eff)}** / {len(have)} gen")
A(f"- `g_obs == 0` esitligi bulunan: "
  f"{sum(1 for r in have if isinstance(r['n_tie_zero'], int) and r['n_tie_zero'] >= 1)}")
A(f"- trivial esitligi (r <= 3 A): "
  f"{sum(1 for r in have if isinstance(r['n_tie_trivial'], int) and r['n_tie_trivial'] >= 1)}")
A(f"- `n_zero >= 1` (g_obs == 0 kabuk): "
  f"{sum(1 for r in have if isinstance(r['n_zero'], int) and r['n_zero'] >= 1)}")
A("")
A("Toplamlar: "
  f"n_tie = {sum(r['n_tie'] for r in have if isinstance(r['n_tie'], int))}, "
  f"bunun g_obs>1 olani = {sum(r['n_tie_g_gt_1'] for r in have if isinstance(r['n_tie_g_gt_1'], int))}, "
  f"g_obs==0 olani = {sum(r['n_tie_zero'] for r in have if isinstance(r['n_tie_zero'], int))}")
A("")
if tie_eff:
    A("### E'yi etkileyebilen esitlikler (g_obs > 1)")
    A("")
    A("| gene | branch | N_P | n_shells | n_tie | n_tie_g_gt_1 | tie_radii_g_gt_1 | domain_changed |")
    A("|---|---|---:|---:|---:|---:|---|:--:|")
    for r in sorted(tie_eff, key=lambda x: -x["n_tie_g_gt_1"]):
        A(f"| `{r['gene']}` | {r['branch_taken']} | {r['N_P']} | {r['n_shells']} | "
          f"{r['n_tie']} | **{r['n_tie_g_gt_1']}** | {r['tie_radii_g_gt_1']} | "
          f"{r['domain_changed']} |")
    A("")
A("### Tum esitlikler, gen basina (r > 3 A)")
A("")
A("| gene | branch | N_P | n_shells | n_tie | g>1 | g==0 | tie_radii |")
A("|---|---|---:|---:|---:|---:|---:|---|")
for r in sorted(tie1, key=lambda x: -x["n_tie"]):
    tr = r["tie_radii"]
    if len(tr) > 60: tr = tr[:60] + f"... (+{tr[60:].count(';')+1})"
    A(f"| `{r['gene']}` | {r['branch_taken']} | {r['N_P']} | {r['n_shells']} | "
      f"{r['n_tie']} | {r['n_tie_g_gt_1']} | {r['n_tie_zero']} | {tr} |")
A("")
A("## Domain degisimi")
A("")
A(f"- `domain_changed = TRUE`, tum genler: **{len(ch)}** / {len(have)}")
A(f"- `domain_changed = TRUE`, yalnizca `pcf_derived` (kuralin isirdigi genler): "
  f"**{len(ch_pcf)}** / {len(pcfd)}")
A(f"- `fallback_envelope` genleri ({len(fb)}): domain zarftan gelir, kural isirmaz")
A("")
A("Bos domain etiketleri (eksik dosya DEGIL): "
  f"`E_BOS` domain_gt = {sum(1 for r in have if r['domain_gt'] == 'E_BOS')}, "
  f"domain_ge = {sum(1 for r in have if r['domain_ge'] == 'E_BOS')}; "
  f"`KIRPMA_BOS` domain_gt = {sum(1 for r in have if r['domain_gt'] == 'KIRPMA_BOS')}, "
  f"domain_ge = {sum(1 for r in have if r['domain_ge'] == 'KIRPMA_BOS')}")
A("")
if ch:
    A("| gene | branch | domain_gt | domain_ge | pivot yaricap(lar) | r_hot | r_hot in domain_ge |")
    A("|---|---|---|---|---|---:|:--:|")
    for r in sorted(ch, key=lambda x: x["gene"]):
        A(f"| `{r['gene']}` | {r['branch_taken']} | {r['domain_gt']} | {r['domain_ge']} | "
          f"{r['pivot_radii'] or '—'} | {r['selected_r_hot']} | {r['r_hot_in_domain_ge']} |")
    A("")
else:
    A("`domain_gt` ile `domain_ge` hicbir gende farkli degil.")
    A("")
pv = [r for r in have if r["pivot_radii"]]
A(f"- `E_ge`'ye giren ama `E_gt`'de olmayan yaricabi olan gen: **{len(pv)}**")
if pv:
    A("")
    A("| gene | branch | pivot yaricap(lar) | domain degisti mi |")
    A("|---|---|---|:--:|")
    for r in sorted(pv, key=lambda x: x["gene"]):
        A(f"| `{r['gene']}` | {r['branch_taken']} | {r['pivot_radii']} | {r['domain_changed']} |")
A("")
A("## r_hot kapsama")
A("")
A("`domain_ge` yalnizca bu analizin turettigi aralik; `fallback_envelope` genlerinde "
  "boru hatti r_hot'u [5,25] zarfindan secer, dolayisiyla FALSE beklenen durumdur.")
A("")
inb = [r for r in have if r["r_hot_in_domain_ge"] == "TRUE"]
outb = [r for r in have if r["r_hot_in_domain_ge"] == "FALSE"]
und = [r for r in have if str(r["r_hot_in_domain_ge"]).startswith("TANIMSIZ")]
nor = [r for r in have if r["r_hot_in_domain_ge"] == MISSING]
A(f"- `domain_ge` icinde: **{len(inb)}**")
A(f"- disinda: **{len(outb)}**")
A(f"  - bunlardan pcf_derived: "
  f"{len([r for r in outb if r['branch_taken'] == 'pcf_derived'])}"
  + (f" — {', '.join('`'+r['gene']+'`' for r in outb if r['branch_taken'] == 'pcf_derived')}"
     if any(r['branch_taken'] == 'pcf_derived' for r in outb) else ""))
A(f"  - fallback_envelope: {len([r for r in outb if r['branch_taken'] == 'fallback_envelope'])}")
A(f"- `domain_ge` tanimsiz (E_BOS / KIRPMA_BOS): **{len(und)}**"
  + (f" — {', '.join('`'+r['gene']+'`' for r in und)}" if und else ""))
A(f"- secilen r_hot yok (handoff_02 yok veya r_hot NA): **{len(nor)}**"
  + (f" — {', '.join('`'+r['gene']+'`' for r in nor)}" if nor else ""))
A("")
A("## N_P ile iliski (Spearman)")
A("")
base = [r for r in have if isinstance(r["N_P"], int)]
for lab, key in (("n_tie", "n_tie"), ("n_zero", "n_zero")):
    s = [r for r in base if isinstance(r[key], int)]
    rho = spearman([r["N_P"] for r in s], [r[key] for r in s])
    A(f"- rho(N_P, {lab}) = " + (f"**{rho:+.3f}**" if rho is not None else "tanimsiz")
      + f"  (n = {len(s)})")
s = [r for r in base if isinstance(r["n_pairs"], int) and isinstance(r["n_zero"], int)]
rho = spearman([r["n_pairs"] for r in s], [r["n_zero"] for r in s])
A(f"- rho(n_pairs, n_zero) = " + (f"**{rho:+.3f}**" if rho is not None else "tanimsiz")
  + f"  (n = {len(s)})")
A("")
q = [r for r in base if isinstance(r["n_zero"], int)]
if q:
    qs = sorted(q, key=lambda r: r["N_P"])
    cut = max(1, len(qs) // 4)
    lo_q, hi_q = qs[:cut], qs[-cut:]
    A(f"- en kucuk {len(lo_q)} kohort (N_P {lo_q[0]['N_P']}-{lo_q[-1]['N_P']}): "
      f"medyan n_zero {st.median([r['n_zero'] for r in lo_q]):g}, "
      f"medyan n_tie {st.median([r['n_tie'] for r in lo_q]):g}")
    A(f"- en buyuk {len(hi_q)} kohort (N_P {hi_q[0]['N_P']}-{hi_q[-1]['N_P']}): "
      f"medyan n_zero {st.median([r['n_zero'] for r in hi_q]):g}, "
      f"medyan n_tie {st.median([r['n_tie'] for r in hi_q]):g}")
A("")
A("## MEFV")
A("")
m = next((r for r in rows if r["gene"] == "MEFV"), None)
if m is None:
    A(f"MEFV `{ROOT}/` altinda yok.")
else:
    A(f"- branch: {m['branch_taken']} · N_P: {m['N_P']} · n_shells: {m['n_shells']}")
    A(f"- n_tie: {m['n_tie']} · tie_radii: {m['tie_radii'] or '—'} · "
      f"n_tie_trivial: {m['n_tie_trivial']} ({m['tie_radii_trivial'] or '—'})")
    A(f"- n_zero: {m['n_zero']} · zero_radii: {m['zero_radii'] or '—'}")
    A(f"- domain_gt {m['domain_gt']} · domain_ge {m['domain_ge']} · "
      f"changed {m['domain_changed']}")
    A(f"- pivot yaricap(lar): {m['pivot_radii'] or '—'}")
    fr11 = genes.get("MEFV")
    if fr11:
        pth = f"{fr11}/04_GLOBAL_CLUSTERING/pair_correlation_plp.tsv"
        if os.path.exists(pth):
            with open(pth) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if fl(r["radius_A"]) is not None and abs(fl(r["radius_A"]) - 11.0) <= EPS:
                        go, p95 = fl(r["g_obs"]), fl(r["g_null_p95"])
                        A(f"- 11 A satiri: g_obs = {go!r} · g_null_p95 = {p95!r} · "
                          f"|fark| = {abs(go - p95):.6g} · in_set_E = {r.get('in_set_E')}")
                        break
    tr = [x for x in (m["tie_radii"] or "").split(";") if x]
    A(f"- beklenen pivot 11 A: 11 esitlik listesinde " +
      ("**VAR**" if "11" in tr else "**YOK**") +
      "; 11 pivot listesinde " +
      ("**VAR**" if "11" in [x for x in (m['pivot_radii'] or '').split(';') if x] else "**YOK**"))
A("")
if nodata:
    A("## VERI_YOK")
    A("")
    A("| gene | eksik |")
    A("|---|---|")
    for r in nodata:
        A(f"| `{r['gene']}` | {r['missing_files'] or 'bilinmiyor'} |")
    A("")
mf = [r for r in have if r["missing_files"]]
if mf:
    A("## Kismi eksik (PCF tablosu var, baska dosya yok)")
    A("")
    A("| gene | eksik |")
    A("|---|---|")
    for r in mf: A(f"| `{r['gene']}` | {r['missing_files']} |")
    A("")

out_md = f"{ROOT}/SHELL_TIE_AUDIT_v3.md"
with open(out_md, "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"  {out_tsv}  ({len(rows)} satir)")
print(f"  {out_md}  ({len(L)} satir)")
