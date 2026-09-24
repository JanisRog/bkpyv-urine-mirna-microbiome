# Final statistical and figure pass after DNA alignment

This is a prospective work plan for the **remaining revision**, written after the first miRNA/MetaPhlAn reanalysis and the 22-sample Kraken2 screen, while the study-wide BKPyV/JCPyV competitive alignment is running. It does not replace the 22 September `ANALYSIS_PLAN.md`, which records decisions made before the initial revised differential results. Any change after seeing the alignment output must be labeled as exploratory.

## Freeze and validate inputs

1. Use the verified 22-recipient accession-to-sample map and blood-DNAemia case/control labels throughout. Preserve the distinct analytical populations: all 22 DNA samples; 16 samples (nine cases, seven controls) for the working human-miRNA model; and the stated sensitivity subsets. The control label does not imply BKPyV-negative urine.
2. Validate all 22 competitive-alignment archives against the fixed BKPyV and JCPyV reference SHA-256, reference lengths, sample IDs, depth files and high-MAPQ read counts. Keep the Kraken2 species assignments and independent genome alignments separate in output tables.
3. Freeze checksums for unique-mature human miRNA counts, viral-miRNA counts, 22 MetaPhlAn profiles, Kraken2 reports, alignment summaries, sample mapping and all analysis code. Retain upstream logs and software/reference version records. Do not invent unavailable extraction blanks, urine volumes, urine PCR results or batch variables.

## Recompute the results from one manifest

- **Human miRNAs:** rerun the recorded DESeq2 main model and edgeR normalization sensitivity from integer unique-mature counts. Show tested feature counts, effect sizes, confidence intervals and adjusted p-values. Repeat all previously specified sample-inclusion and S3/S24 influence analyses without choosing a subset by significance. No claim based on nominal p-values alone.
- **Viral miRNAs:** report raw and two normalized views, explicitly naming the 3p feature as shared BKPyV/JCPyV. Compare with the DNA genome-alignment evidence by sample without treating RNA and DNA library counts as directly commensurate or an absolute viral-load measurement.
- **Bacteria:** rebuild class, family, genus and species tables from all 22 MetaPhlAn profiles. Make the all-22 DNA result primary and use the paired miRNA subsets only as sensitivities. Recompute compositional comparisons, FDR and full numeric tables; report the range and individual sample values that a group mean conceals. Do not equate bacterial-only relative percentages with fractions of all DNA reads.
- **Viral DNA:** use Kraken2 for broad screening and high-MAPQ competitive alignment plus genome breadth for BKPyV/JCPyV species support. Describe group distributions and exceptions such as S24 using the correct blood-DNAemia labels. Treat reference-mapping evidence as DNA detection, not productive infection or copies/mL.
- **Fungi, protozoa and residual human DNA:** supply per-sample Kraken2 counts, fractions of classifier input and distinct minimizers. These are screens; do not label fungal/protozoan infections or claim an exhaustive parasite assay. Describe residual-human assignments as a post-host-filter classifier estimate, not a complete human DNA percentage. Without available negative controls or quantitative calibration, avoid absolute abundance claims.
- **Integration:** rerun only the fixed historical miRNA–bacterial candidate set, conditioning on clinical group and applying the recorded multiplicity correction. Label any new cross-domain observation prompted by the genome results exploratory.

## Figures and supplementary material

1. **Study and QC figure:** show the 22-sample provenance, the 16-sample human-miRNA subset, human assigned counts, host-filtered DNA pairs and residual-human classifier fractions. Annotate unusually low-yield samples without silently excluding them.
2. **Human-miRNA result figure:** keep an effect-size plot for the historical candidates with confidence intervals, numeric effect/FDR values in a companion table, and clear labeling of untested features. Do not imply a discovery where none survives FDR.
3. **Virome figure:** paired per-sample BKPyV/JCPyV DNA evidence (high-MAPQ first-mate counts and genome breadth) alongside separately scaled viral-miRNA counts. Mark cases by blood DNAemia and make S14/S24 legible. Do not merge RNA and DNA scales.
4. **Bacterial figure:** show individual-level class distributions and a readable family/genus/species view of the taxa discussed clinically. Prefer dot/interval or compact heatmap views over using a crowded 22-bar stack as the only evidence. Give exact group summaries and adjusted tests in a table.
5. **Sankey:** retain the user-preferred descriptive Sankey. Check every displayed class has visible links, keep left and right scales explicit, label the shared BK/JCPyV 3p feature, and state that links indicate sample identity rather than biological interaction. Put a static version in the manuscript and the interactive version in the author package if journal policy permits.
6. **Supplement:** full tested human-miRNA table, four-rank bacterial tables, sample-level viral DNA and RNA counts, fungal/protozoan screen, DNA/miRNA QC, sensitivity analyses, accession mapping and software/reference provenance.

## Manuscript and response pass

Regenerate every number and plot from the frozen inputs, then update the Results, expanded Discussion, figure legends and point-by-point response together. Specifically explain the distinction between blood DNAemia and urinary shedding, the shared 3p miRNA, discordance between original and revised analyses, low-biomass contamination limits, and why absolute reads or copies per mL cannot be reconstructed. Audit each claim against its table and figure. Render and visually check the final DOCX files. Keep the repository private until the author approves a release.
