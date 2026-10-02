# Joint human and viral miRNA count analysis

The primary comparison is nine blood-BKPyV-DNAemia cases versus seven controls with at least 2,000 assigned human mature-miRNA counts. It contains 244 filtered human features and two prespecified viral sequences. DESeq2 fits all 246 count rows in one negative-binomial group model. Size factors come from the earlier human-only positive-count calculation and are fixed before the viral rows enter the model; thus high viral counts cannot change the exposure denominator. Dispersion is fit locally, automatic outlier replacement and independent filtering are off, and Cook's-distance filtering is on. Benjamini–Hochberg correction covers every finite human and viral test in each fit. Robust edgeR quasi-likelihood models reuse the host-only library totals and TMMwsp factors. This is the common statistical basis for Figure 1; raw- and retained-read-scaled permutations remain denominator sensitivity checks, not the primary tests.

| Primary 9/7 result | DESeq2 log2FC | DESeq2 p | Joint q | Robust edgeR q |
| --- | ---: | ---: | ---: | ---: |
| BKPyV/JCPyV-shared 3p | 3.495 | 0.000344 | 0.057 | 0.198 |
| BKPyV 5p | 4.758 | unavailable | unavailable | 0.198 |

BKPyV 5p is not placed at a numerical p-value ordinate in Figure 1. S7 has a normalized 5p count of about 53,601 and Cook's distance 7.127; the DESeq2 primary test is filtered. Turning off Cook's filtering produces a small p-value, but bypasses an outlier safeguard and is diagnostic only. Re-estimating the main model without S7 reduces the 5p log2FC to 2.951 and yields q=0.140. S7 is not removed from the primary analysis.

In the all-22 sensitivity, BKPyV 5p has DESeq2 q=0.012 but robust edgeR q=0.374. The 3p q is 0.065 with DESeq2 and 0.378 with edgeR. Inclusion subsets and statistical method therefore affect viral inference. The 3p mature sequence is identical between BKPyV and JCPyV and cannot be attributed to one species from these small-RNA counts. Every urine DNA library has BKPyV coverage, including controls defined by negative *blood* DNAemia.

The 15 human features meeting the source paper's descriptive display rule (unadjusted p<0.05 and absolute log2FC at least 1) are unchanged by adding the viral rows. Their p and q values change slightly because the joint model re-estimates the dispersion trend and corrects across 231 finite tests; no human feature has q<0.05 in any of the four inclusion strategies. The source report's 8/8 subset and CLC workflow remain different, so this is not an exact replication.

Run `analysis/run_joint_mirna_models.R` after `analysis/run_count_models.R` has produced the host-only factors. `analysis/make_joint_mirna_figure.py` reads its primary result table and writes the figure, plot-data CSV, caption and summary JSON. All source counts and numerical outputs remain in the private local author package, not this code-only repository.
