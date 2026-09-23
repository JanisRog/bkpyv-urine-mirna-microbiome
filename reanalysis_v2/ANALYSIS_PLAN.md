# Revision analysis specification

Frozen before viewing the new differential-expression results, 22 September 2026. This is a retrospective revision, not a preregistration. Original outputs are retained separately.

## Human miRNAs

Input: cluster-verified unique-mature integer counts, preserving read multiplicity while excluding assignments ambiguous across distinct mature names. Original precursor-summed counts are historical comparison only. Fractional assignments remain a sensitivity output and will not be rounded into DESeq2 input.

The working main subset uses at least 2,000 human unique-mature assigned counts (9 cases, 7 controls), without the arbitrary 30-feature requirement. This is a pragmatic, unvalidated yield criterion. Mandatory sensitivity subsets: all 22 samples (11/11), at least 1,000 counts (9/8), and the previous combined 2,000-count/30-feature rule (8/7). No subset will be selected based on significance.

DESeq2 negative-binomial Wald tests, group-only design, positive-count size factors, no automatic count replacement, Cook's outlier filtering retained, and BH correction over finite tested p-values within each analysis. Explicit feature prefilter: at least 10 counts in at least 3 included samples. Independent filtering disabled to make the tested family explicit. Local dispersion fit specified in advance for these sparse data. Unshrunk log2 fold change, standard error and Wald intervals reported; no significance claim based on nominal p alone. Normalization sensitivity: edgeR TMMwsp with quasi-likelihood testing. Missing batch/clinical data preclude adjusted causal interpretation.

S3 and the highly concentrated S24 will receive explicit leave-one-out influence analyses of the working main subset. These are sensitivities, not data-driven exclusions. Inspect size factors, library-size relationships, expression distributions and sample composition. Methods cannot establish absolute RNA concentrations without an external anchor.

## Viral miRNAs

Preserve validated historical BK-only counts and shared BK/JC 3p attribution caveat. Report per-sample raw counts and sensitivity after scaling by human model size factors and by retained small-RNA reads. Explore case-control differences by label permutations on log2-scaled values; treat this as descriptive/exploratory, with BH over both viral features per scaling/subset. Human normalization assumes comparable host-miRNA background and does not measure absolute viral burden. No negative-binomial dispersion estimate from two viral features alone.

## Bacteria

Use all 22 DNA profiles as the main available table analysis; miRNA QC is not DNA QC. Include each of the three reduced paired subsets as sensitivities. Retain class/family/genus/species summaries and explicit compositional denominators. Relative percentages are not read counts; do not supply them to count-based ALDEx2/DESeq2. Exploratory group tests will assess CLR-transformed relative profiles with explicit zero replacement sensitivity, alongside relative-percentage contrasts. BH within each rank and globally across rank/taxon tests. No contamination correction or absolute abundance is claimed from these tables. Upstream DNA QC and broad-domain classifications remain pending cluster records or reruns.

## Integration and presentation

Historical 12 candidates are fixed for comparison, not chosen anew from nominal significance. Associations with bacterial classes are exploratory; group-conditioned analysis and FDR must address shared group structure and multiplicity. Sankey remains descriptive, with explicit normalization, separate domain scales and display-only grouping. All plots, tables and manuscript claims must identify the sample subset and counting/normalization policy.

## Unresolved requirements

Exact original miRBase/miRDeep2 versions and input-read history; clinical covariates and the symptomatic E. coli control's ID; urine fraction and input volume; DNA trimming, host-filter/QC and database provenance; viral/fungal/protist read-level analysis. No reanalysis can manufacture missing blanks, biological replicates or quantitative calibration. These gaps must be disclosed rather than represented as completed.
