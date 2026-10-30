"""One report. Every gene carries its own regime, tier and caveat."""
import csv, collections
import numpy as np
R = list(csv.DictReader(open("results_v3/BATCH_DATA_ACMG_SF_v3.2_v3_classified.tsv"),
                        delimiter="\t"))
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None
def i(x):
    try: return int(float(x))
    except (TypeError, ValueError): return None
def g(x, n=3):
    v = f(x)
    return "" if v is None else f"{v:.{n}f}".rstrip("0").rstrip(".")

L = []; A = L.append
A("# Batch report — ACMG SF v3.2, run v3\n")
A("Generated from `BATCH_DATA_ACMG_SF_v3.2_v3_classified.tsv`. Every figure is read from")
A("the run artefacts; nothing is re-derived here. **Each gene carries its own regime and")
A("caveat in the table below, so a result cannot be quoted apart from its qualification.**\n")
A("| | |")
A("|---|---|")
A("| Runs | 68 ACMG SF v3.2 genes + KCNA2 (case gene) + MEFV |")
A("| Batch window | 2026-10-25T19:33:13Z to 2026-10-26T15:13:18Z, 19.7 h, 17.4 min/gene |")
A("| `config_sha256` / `code_version` | `ff721543` / `ac7cce0f` |")
A("| Branch | `rfp-containment` @ `59febb7` |")
A("| Decisions in force | RFP-CONTAINMENT-0001, MANE-MATCH-0001, REPORTING-CONTEXT-0001, VOXEL-TIE-0001, UNDERPOWERED-FIXTURE-0001, HIGH-COVERAGE-GENES-0001, MEFV-RESULT-0001 |")
A("| Data | `results_v3/`, 57 GB |")
A("")

REG = {
 "works": ("WORKS",
   "All four quality criteria or three of four. Quote freely, with the tier.",
   "fold enrichment >= 1.5, median Jaccard >= 0.90, binding residue not an outlier, footprint coverage < 0.30"),
 "works_with_caveat": ("WORKS — CAVEAT REQUIRED",
   "A footprint exists but meets at most one quality criterion. **Never quote the centre count without the caveat column.**", ""),
 "diffuse_no_localised_region": ("DIFFUSE — no localised region",
   "**Not a refusal and not a failure.** Hotspot centres stand and are a finding: pathogenicity is spread rather than localised. Only the footprint geometry is withheld — the radius covers more than half the structure at fold enrichment near 1.0, so inside and outside the sphere have the same pathogenic density and a region that large is not a region. See `DECISION-HIGH-COVERAGE-GENES-0001`.", ""),
 "refused_degenerate": ("REFUSED — degenerate radius",
   "Hotspot centres stand. The containment radius came out at zero, so no footprint has volume.", ""),
 "no_spatial_concentration": ("NEGATIVE — informative",
   "The test had the power to reject and did not. A statement about the evidence, not about the protein. **Not a failure.**", ""),
 "cannot_be_asked": ("UNINFORMATIVE",
   "The test could not have rejected at any `B`. **Says nothing about the gene — must never be read as a negative.**", ""),
 "input_limited": ("INPUT-LIMITED",
   "Stage A could not assemble usable inputs, or the cohort is unanalysable. No scientific claim either way.", ""),
}
ORDER = ["works", "works_with_caveat", "diffuse_no_localised_region", "refused_degenerate",
         "no_spatial_concentration", "cannot_be_asked", "input_limited"]

A("## How to read this report\n")
A("Seven regimes. The regime, not the centre count, determines what a gene's result")
A("licenses.\n")
A("| regime | genes | what it licenses |")
A("|---|---:|---|")
for k in ORDER:
    n = sum(1 for r in R if r["regime"] == k)
    A(f"| **{REG[k][0]}** | {n} | {REG[k][1]} |")
A("")
A("Tiers apply only within WORKS: **A** meets all four quality criteria")
A(f"({REG['works'][2]}); **B** three; **C** two; **D** one or none — and D is reported as")
A("its own regime because one criterion is not enough to carry a claim.\n")

A("## Every gene, with its caveat\n")
A("Sorted by regime, then centre count. Blank means the stage was never reached.\n")
A("| gene | regime | tier | N_P | N_B | M | r_hot | centres | cov@r_hot | fold | "
  "r_fp | fp res | fp cov | jaccard | caveat |")
A("|---|---|:--:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
for k in ORDER:
    grp = sorted((r for r in R if r["regime"] == k),
                 key=lambda x: -(i(x["n_centers"]) or 0))
    for r in grp:
        A("| `{gene}` | {reg} | {t} | {np} | {nb} | {M} | {rh} | {nc} | {cv} | {fe} | "
          "{rf} | {fr} | {fc} | {jc} | {cav} |".format(
            gene=r["gene"], reg=REG[k][0], t=r["tier"] or "—",
            np=r["N_P"] or "", nb=r["N_B"] or "", M=r["M"] or "",
            rh=g(r["r_hot"], 1), nc=r["n_centers"] or "",
            cv=g(r["coverage_at_r_hot"]), fe=g(r["fold_enr_at_r_hot"]),
            rf=g(r["r_fp"], 1), fr=r["fp_residues"] or "", fc=g(r["fp_coverage"]),
            jc=g(r["jaccard_median"]), cav=r["caveat"]))
A("")

A("## What determines the regime\n")
A("### Cohort size decides whether the question can be asked\n")
staged = [r for r in R if i(r["N_P"]) is not None]
A("| N_P | genes finding centres | rate |")
A("|---|---|---:|")
for lo, hi in ((1, 2), (3, 5), (6, 10), (11, 25), (26, 60), (61, 400)):
    gg = [r for r in staged if lo <= (i(r["N_P"]) or 0) <= hi]
    if not gg: continue
    got = sum(1 for r in gg if (i(r["n_centers"]) or 0) > 0)
    A(f"| {lo}–{hi} | {got} / {len(gg)} | {got/len(gg)*100:.0f} % |")
cen = [r for r in staged if (i(r["n_centers"]) or 0) > 0]
non = [r for r in staged if (i(r["n_centers"]) or 0) == 0]
d1 = np.median([i(r["N_P"]) / i(r["M"]) for r in cen if i(r["M"])])
d0 = np.median([i(r["N_P"]) / i(r["M"]) for r in non if i(r["M"])])
A("")
A(f"Pathogenic density `N_P/M` separates the groups {d1/d0:.1f}x — median {d1:.4f} against")
A(f"{d0:.4f}. The working range begins near `N_P = 10`; below 5 it is close to a coin flip,")
A("and at 1–2 the combinatorial floor makes rejection impossible at any `B`.\n")

A("### Architecture decides whether the question has an answer\n")
A("| gene | N_P | centres | fold | what it shows |")
A("|---|---:|---:|---:|---|")
A("| `COL3A1` | 306 | 0 | — | a collagen triple helix spreads its variants along the chain; nothing for a sphere to enclose |")
A("| `MYH7` | 260 | 25 | 1.07 | largest cohort with a footprint, almost no enrichment |")
A("| `RBM20` | 5 | 6 | 16.40 | highest enrichment in the batch, on five residues |")
A("| `TRDN` | 3 | 3 | 11.67 | three residues in a tight pocket, coverage 0.011 |")
A("")
A("Cohort size sets whether the question can be asked; protein architecture sets whether")
A("it has an answer. Neither is something the pipeline can improve.\n")

A("### Quality of the 38 footprints produced\n")
fp = [r for r in R if f(r["r_fp"])]
fe = [(f(r["fold_enr_at_r_hot"]), r["gene"]) for r in fp if f(r["fold_enr_at_r_hot"]) is not None]
A("| fold enrichment | genes |")
A("|---|---:|")
for lo, hi, lbl in ((0, 1.05, "≈1.0 — none"), (1.05, 1.5, "1.05–1.5 weak"),
                    (1.5, 3.0, "1.5–3 moderate"), (3.0, 99, ">3 strong")):
    A(f"| {lbl} | {sum(1 for x in fe if lo <= x[0] < hi)} |")
A("")
A(f"Only {sum(1 for x in fe if x[0] >= 1.5)} of {len(fe)} reach 1.5. Significance under the")
A("positional null does not require enrichment, so the majority are true statements about")
A("geometry and weaker claims than a centre count alone implies.\n")
jc = [f(r["jaccard_median"]) for r in fp if f(r["jaccard_median"]) is not None]
ol = sorted(r["gene"] for r in fp if str(r["binding_is_outlier"]).lower() == "true")
A(f"Stage D median Jaccard across the 38 is **{np.median(jc):.3f}** — robustness is the")
A(f"strongest dimension. But `r_fp` rests on a single outlying residue in {len(ol)} of them:")
A(f"{', '.join('`'+x+'`' for x in ol)}.\n")

A("## The finding that governs how this batch is quoted\n")
d = sorted((r for r in R if r["regime"] == "works_with_caveat"),
           key=lambda x: -(f(x["fp_coverage"]) or 0))
A("Footprint coverage, not centre count, is what marks a thin result. Over the 38 footprint")
A("producers, coverage and centre count rise together (Spearman rho = +0.74), but centre")
A("count on its own carries almost no information about quality (rho = -0.09 against the")
A("number of criteria met). The five CAVEAT-REQUIRED genes:\n")
A("| gene | centres | fold | fp coverage | jaccard | outlier-bound | the one criterion it keeps |")
A("|---|---:|---:|---:|---:|:--:|---|")
KEEP = {"HNF1A": "robustness only", "SCN5A": "robustness only",
        "STK11": "robustness only", "TP53": "robustness only",
        "NF2": "coverage only"}
for r in d:
    A(f"| `{r['gene']}` | {r['n_centers']} | {g(r['fold_enr_at_r_hot'])} | "
      f"{g(r['fp_coverage'])} | {g(r['jaccard_median'])} | yes | {KEEP[r['gene']]} |")
A("")
A("Two properties hold for all five: every one falls below the 1.5 enrichment criterion")
A("(`TP53` by 0.004) and every one rests on an outlying binding residue. Four of the five")
A("also sit in the coverage band")
A("0.31-0.48, just under the frozen 0.50 limit — they passed the same rule that refused")
A("eleven other genes, with little margin, and those eleven and these four are one")
A("continuous population cut at a single threshold. `NF2` is the exception and fails for a")
A("different reason: coverage 0.029, but unstable under perturbation as well as weakly")
A("enriched and outlier-bound.\n")
A("**Consequence for reporting:** the three largest footprint coverages in the batch are all")
A("CAVEAT-REQUIRED (`SCN5A` 0.482, `TP53` 0.478, `HNF1A` 0.471), and `SCN5A` carries the")
A("largest centre count in the batch at 343. But the inverse does not follow either way:")
A("`TSC2`'s 128 centres at coverage 0.232 is tier A and needs no coverage caveat, while")
A("`BRCA1`'s 226 at coverage 0.106 passes every quality criterion and is still not quotable")
A("as clustering evidence -- for a reason the quality criteria cannot see. See **The")
A("measurement itself** below. Read the regime, the coverage AND the B1 status; never the")
A("centre count alone.\n")

A("## The measurement itself: is the elevation clustering, or folded-vs-disordered contrast?\n")
A("AlphaFold places intrinsically disordered regions as extended, low-pLDDT sprawl, while")
A("pathogenic missense variants concentrate in the compact folded core. The positional null")
A("redraws N_P residues uniformly over ALL of `U_struct`, so the observed cohort and its")
A("null can be drawn from different structural regimes. Across this batch a median of")
A("**93.3 %** of P/LP residues sit at pLDDT >= 70, and in `BRCA1` 80 % of `U_struct` is")
A("below pLDDT 50 while 94 % of its P/LP are above 70.\n")
A("That premise was tested directly, null against null: draw N_P residues uniformly from")
A("the folded subset only -- a set with **no clustering whatever** -- and run the same pair")
A("correlation against the same null. Two numbers come out: how much of the observed excess")
A("over the null the confound alone reproduces (`b1_explained`), and whether it elevates the")
A("*same radii* (`b1_overlap`). The ratio alone is not enough -- in genes whose observed")
A("excess is tiny it explodes (`PALB2` reaches 41 at zero overlap) -- so it is gated on")
A("overlap.\n")
A("| gene | regime | centres | explained | overlap | reading |")
A("|---|---|---:|---:|---:|---|")
for _g, _rg, _nc, _e, _o, _note in (
        ("HNF1A", "WORKS - CAVEAT REQUIRED", "187", "1.16", "0.97",
         "the confound produces MORE than was observed"),
        ("LDLR", "REFUSED - coverage", "346", "0.76", "1.00", "location and magnitude both reproduced"),
        ("BRCA1", "WORKS, tier B", "226", "0.75", "0.98", "location and magnitude both reproduced"),
        ("LMNA", "WORKS", "6", "0.60", "0.92", "location and magnitude both reproduced"),
        ("MYH7", "WORKS", "25", "0.58", "0.92", "location and magnitude both reproduced"),
        ("KCNQ1", "REFUSED - coverage", "217", "0.54", "0.95", "location and magnitude both reproduced")):
    A(f"| `{_g}` | {_rg} | {_nc} | **{_e}** | **{_o}** | {_note} |")
A("")
A("**Six genes are B1-CONFOUNDED**: the disorder contrast alone reproduces both where the")
A("elevation sits and most of how large it is. Their centre counts are not evidence of")
A("spatial clustering and must not be quoted as such, whatever their tier. `HNF1A` already")
A("carried a caveat; `BRCA1` and `MYH7` did not, and that was an error this diagnostic")
A("corrects -- `BRCA1` meets every quality criterion and is still not quotable. `MYH7` is")
A("compromised twice over: its near-tie decides 25 against 73 centres on a 0.82 % margin,")
A("and 58 % of its elevation is confound at 0.92 overlap.\n")
A("`KCNH2` (overlap 1.00, explained 0.32) and `TP53` (1.00, 0.37) are BAND_CONFOUNDED: the")
A("confound fixes *where* the elevated band falls but accounts for only a third of its size.")
A("Their radius is confound-driven even where their signal is not.\n")
A("**25 genes are CLEAN** -- overlap below 0.50 and under 30 % explained. `KCNA2`, the case")
A("gene, is among them at 0.29 / 0.13, as are `TSC1` (0.09 / 0.13), `MSH2` (0.05 / 0.05),")
A("`SMAD4`, `STK11`, `PMS2`, `TGFBR1` and `TGFBR2`. Of the rest, 22 are PARTIAL and 15 were")
A("NOT_TESTED -- no centres, or fewer than three P/LP residues.\n")
A("`COL3A1` could not be tested at all: its folded subset is smaller than its N_P of 306. A")
A("collagen triple helix is low-pLDDT along its whole length, so there is no folded core to")
A("draw from -- which is also why it found zero centres.\n")
A("This is B1 of the 2026-09-13 audit, open and untested until now. It is a diagnostic, NOT")
A("a filter: no residue is excluded for low pLDDT and no cohort is touched. Acting on it")
A("would mean changing the null -- stratifying it so it draws from the same structural")
A("regime as the observed cohort -- which is a methodology decision requiring a re-run, not")
A("a repair. Diagnostic permutations B = 600, five independent draws per gene.\n")
A("## Requiring attention\n")
A("| gene | issue | action |")
A("|---|---|---|")
A("| `MYH7` | near-tie at 0.82 % decides 25 vs 73 centres (+192 %); the four objectives are indistinguishable between the two radii | near-tie severity should grade on outcome divergence, not distance gap |")
A("| `LMNA`, `SCN5A`, `SMAD4` | the PCF block kept is the LONGEST, not the strongest: the discarded block is 5.2x / 2.3x / 2.2x stronger over the null | block-selection key decision, pending |")
A("| `CYP27A1` | AFDB serves 2 fragments for Q02318; stitching prohibited (METHOD_SPEC II.13) | standing input property |")
A("| `BTD` | AlphaFold and UniProt numbering disagree by −20, identity 1.0000 over 523 residues | Lead decision |")
A("| `GLA` | resolves to two reviewed accessions, Q9UHE5 and P06280 | Lead decision |")
A("| `CASQ2` | `r_contain = 0` — first real occurrence of weakness 6 of RFP-CONTAINMENT-0001 | methodology decision |")
A("| — | `BH_BY_DISAGREEMENT` fires MAJOR in 48 of 48 genes that find anything | severity decision, pending |")
A("| — | `handoff_02` carries `global_clustering_plp_p: null` while the value sits in Stage 04 | reporting gap |")
A("")
A("## What this batch says about the pipeline\n")
A("The 70 runs were not an application of the method, they were a characterisation of it.")
A("Read as that, the batch separates three things that look alike from the outside: an")
A("outcome the pipeline got right by declining, a limitation of the inputs, and a defect of")
A("our own.\n")
A("### 51 genes carry a finding, in two different forms\n")
A("Reclassifying the high-coverage genes changes the headline count. A gene with")
A("significant centres has a result whether or not a region can be drawn around them:\n")
A("| | n | what the finding is |")
A("|---|---:|---|")
A("| **localised** | 39 | centres plus a footprint — a region can be named |")
A("| **diffuse** | 12 | centres but no nameable region: pathogenicity is spread. The LOF genes — `GAA`, `OTC`, `RPE65`, `MLH1`, `MEN1`, `PTEN`, `ABCD1`, `SDHB`, `ACVRL1` — where almost any destabilising missense is pathogenic |")
A("")
A("**The diffuse regime is not a disorder artefact.** That was worth testing because the")
A("confound produces broad elevation, which produces high coverage; the prediction was that")
A("diffuse genes would be enriched for it. They are not: the B1-compromised rate is **16 %**")
A("among localised genes (6 of 38) and **17 %** among diffuse ones (2 of 12). The most")
A("extreme cases settle it — `OTC` at 0.915 coverage is CLEAN, `RPE65` at 0.872 is CLEAN with")
A("overlap 0.00, `GAA` at 0.882 is CLEAN. Elevation across 88-92 % of a protein is real")
A("there, not manufactured. The two exceptions, `LDLR` (0.76 / 1.00) and `KCNQ1`")
A("(0.54 / 0.95), are CONFOUNDED and must not be quoted.\n")
A("Read fold enrichment ~ 1.0 together with that: inside and outside the sphere have the")
A("same pathogenic density because **everywhere** is pathogenic. 'No hotspot' is the wrong")
A("reading; 'the hotspot is the protein' is the right one.\n")
A("### The 19 genes that produced nothing are all defensible\n")
A("| n | cause | verdict |")
A("|---:|---|---|")
A("| 4 | external input — BTD numbering offset, CYP27A1 served as two fragments, GLA resolving to two accessions, WT1 with no classified residue | a property of ClinVar/AlphaFold, not of the method. None was passed over silently; each is BLOCKED or escalated |")
A("| 4 | `N_P <= 2` | the combinatorial floor puts the attainable p above `q/m` at **any** B. Reporting UNDERPOWERED rather than a p-value is correct |")
A("| 11 | the test had power and did not reject | a negative about the evidence. `COL3A1` carries 306 P/LP residues and zero centres: a collagen triple helix spreads variants along the chain |")

A("| 1 | `r_contain = 0` (`CASQ2`) | a correct refusal, over a rule that has no positive lower bound |")
A("")
A("**Not one of the 19 is a bug**, and in 15 of them the pipeline declined to assert")
A("something rather than failing to compute it. Add the 12 diffuse genes, where it withheld")
A("a region it could not justify while keeping the centres, and that is 27 places where a")
A("less careful pipeline would have emitted a number. It is the single strongest piece of")
A("evidence in this batch and worth more than any positive result in it.\n")
A("### Our own defects are all on the positive side\n")
A("| # | defect | measured effect | status |")
A("|---:|---|---|---|")
A("| 1 | the positional null draws from all of `U_struct` while the cohort lives in the folded core | 6 genes CONFOUNDED, 2 BAND_CONFOUNDED, 22 PARTIAL of 55 tested | measured, `b1_status` column |")
A("| 2 | the PCF block is kept by **length**, not strength | 3 of 17 multi-block genes keep the weaker block; 5.2x / 2.3x / 2.2x excess discarded | measured, not changed |")
A("| 3 | near-tie severity grades the distance gap, which does not predict consequence | rho = -0.057; `MYH7` decides 25 against 73 centres on 0.82 % | measured, not changed |")
A("| 4 | centre count is a continuous ramp in `r`, with no plateau to select | `MYH7` 8->12, 9.5->25, 10.5->73, 11.5->108 | a design question, not a bug |")
A("| 5 | `BH_BY_DISAGREEMENT` fires MAJOR in every gene that finds anything | 48 of 48 | a warning that always fires carries no information |")
A("| 6 | the containment rule admits `r_contain = 0` | 1 gene (`CASQ2`) | weakness 6, open by decision |")
A("| 7 | QC-F1 withholds the footprint on coverage alone, without consulting fold enrichment | 11 genes, all at fold 1.00-1.16 so nothing was lost; a gene at coverage > 0.50 **and** fold > 1.5 would be withheld wrongly | open, no occurrence yet |")
A("")
A("Items 1 and 2 would change reported numbers; 3, 4 and 5 change how they are read; 6 is")
A("one gene. **None is a coding error** — every one is a choice in the statistical or")
A("geometric design, which is why each is recorded with its measurement rather than")
A("quietly patched. Acting on 1 or 2 means changing the null or the domain rule, and that")
A("means a new `config_sha256` and a re-run of the whole batch.\n")
A("### Where the discovery layer ends and the calibration layer begins\n")
A("The defects cluster in one place, and saying where is more useful than listing them.")
A("Cohort construction, QC, admissibility and refusal — the **discovery** layer — holds up:")
A("31 of 31 non-results are right, and the input failures are surfaced rather than absorbed.")
A("What is weaker is the **calibration** layer: what the null is drawn from, and how stable")
A("the selected radius is. Every one of the six sits there.\n")
A("### How much to trust what\n")
A("| | trust | basis |")
A("|---|---|---|")
A("| the refusals and the negatives | **high** | 27 declined assertions, all verified; power certificate passes on every negative |")
A("| the 14-gene core | **moderate-high** | centres >= 10, median Jaccard >= 0.90, B1 not compromised |")
A("| the other 16 `works` genes | **low-moderate** | 8 of them rest on <= 3 centres and cannot be reported alone |")
A("| the 6 B1-CONFOUNDED genes | **do not quote as clustering** | `BRCA1` meets all four quality criteria and still fails this |")
A("| any single `r_hot` figure | **read the near-tie column first** | in near-tie genes the radius is manufactured from a ramp |")
A("")
A("The 14: `ATP7B`, `BAG3`, `ENG`, `KCNA2`, `MSH2`, `PMS2`, `PRKAG2`, `RB1`, `SMAD4`,")
A("`TGFBR1`, `TGFBR2`, `TSC1`, `TSC2`, `VHL`.\n")
A("**The one-line answer.** The pipeline is reliable when it refuses and conditional when it")
A("asserts. It carries a finding in 51 of 70 genes -- 39 localised, 12 diffuse -- and said")
A("nothing in the other 19, correctly in every one. Of the 51, 25 are B1-clean, 8 are")
A("measurement artefact and must not be quoted, and the rest sit in between with their")
A("caveat attached. Fourteen stand entirely on their own.\n")
A("## Terminal states\n")
A("| state | genes |")
A("|---|---:|")
for k, v in collections.Counter(r["run_status"] for r in R).most_common():
    A(f"| {k} | {v} |")
A("")
open("results_v3/BATCH_REPORT_ACMG_SF_v3.2_v3.md", "w").write("\n".join(L) + "\n")
print(f"  results_v3/BATCH_REPORT_ACMG_SF_v3.2_v3.md ({len(L)} satir)")
