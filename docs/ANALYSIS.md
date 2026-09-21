# Analysis scope and outputs

## Reproduced analysis

Samples require at least 2,000 total human-miRNA counts and at least 30 human features with 10 or more counts. Human features require at least 10 counts in at least three retained samples. Viral counts do not enter sample QC.

Human counts are transformed as log2(count + 1), then quantile-normalized with the original minimum-rank tie assignment. Viral counts are transformed as log2(count + 1) without size normalization. Two-sided Welch tests compare transformed group means. Candidate thresholds are nominal p ≤ 0.05 and absolute transformed mean difference ≥ 1. Benjamini–Hochberg FDR covers all 268 original tests. Transformed mean differences are not measured fold changes in copies per mL.

Human normalization sensitivity averages reference quantiles across the ranks occupied by ties. Viral sensitivity scales each count to one million total human-miRNA counts before log2(value + 1). These sensitivities isolate particular analytical choices; neither establishes an optimal final primary model.

Bacterial profiles are closed to 100% and aggregated to class, family, genus and species. Group averages weight each sample equally. Exact permutation tests enumerate the 11,440 allocations of nine case labels to 16 profiles and compare absolute mean percentage differences. These observational tests assume exchangeability and do not adjust for clinical covariates or compositional dependence. FDR is reported within rank and across all 184 rank–taxon tests. All-22 comparisons are descriptive.

The reproduced correlations use raw-miRNA-count ranks and bacterial percentages across the selected samples. Sensitivities use human quantile-normalized values and human-count-scaled viral values. FDR covers the full correlation matrix. Within-group correlations are descriptive only. Candidate selection on the same cohort and shared case/control structure limit interpretation.

## Sankey

`Figure_3_sankey` retains the original numerical link definitions:

- Left: class relative percentage in each sample. Each sample receives 100 percentage points in total; summed class node heights are not group mean percentages.
- Right: absolute difference between a sample's log2(raw miRNA count + 1) and the mean log2 value of the seven controls. Only absolute differences ≥ 0.5 are drawn. Each control is included in its reference mean.
- Each side uses one linear width conversion shared by all its links, with separate width keys and a visible break. There is no sample-specific rescaling or conservation requirement between sides.
- Right colors indicate sample-specific direction: red/green for human above/below the control mean and dark/light blue for viral above/below. Neutral miRNA nodes make no significance claim.

The browser HTML is self-contained and exposes exact link values on hover. The 0.5 cutoff is for display; omitted links need not indicate absent expression. The figure is sensitive to raw count yield and does not measure biological interactions.

## Generated files

All outputs remain local and are ignored by Git.

| Output | Contents |
| --- | --- |
| `sample_qc.csv`, `human_miRNA_depth.csv` | Per-sample filters, retained status and count yield. |
| `miRNA_all_tests_reproduced.csv`, `miRNA_candidates.csv` | Original tests, transformed differences, global FDR and candidate flags. |
| `human_average_ties_sensitivity.csv` | Same human features with averaged ties and human-feature FDR. |
| `viral_per_million_human_miRNA_sensitivity.csv` | Human-count-denominator sensitivity; its two-test FDR is distinct from the original global correction. |
| `*_per_sample_16.csv`, `*_summary_16.csv` | Individual bacterial percentages and group mean, median, quartiles, nonzero detection counts and within-rank tests. |
| `*_summary_all22.csv`, `taxa_tests_all_ranks.csv` | Full-cohort descriptive summaries and combined tests with across-rank FDR. |
| `correlations_original_scale.csv`, `correlations_normalized_sensitivity.csv` | Original and normalization-sensitivity rho, p and adjusted p values. |
| `correlations_within_*.csv` | Within-group rho; constant pairs remain missing. |
| `original_paper_candidate_comparison.csv` | Candidate-name/direction comparison with the original abstract; not a reconstruction of its raw analysis. |
| `sankey_links.csv`, `sankey_signed_log2_deviations.csv` | Every plotted link with its units, and signed deviations before display filtering. |
| `run_summary.json` | Numerical summary, versions, input hashes and reproduction status. |

Main figures are `Figure_1_miRNA`, `Figure_2_bacteria` and `Figure_3_sankey`. Supplementary Figures S1–S3 are `Supplement_Family`, `Supplement_Genus` and `Supplement_Species`; S4 is `Supplement_Correlations`. PNG and editable-text SVG are exported, with HTML additionally available for the Sankey.
