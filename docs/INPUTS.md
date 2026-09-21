# Required input tables

Supply an approved local data directory with the eight relative paths below. CSV files use a header row and comma delimiters; UTF-8 with or without BOM is accepted by the pandas loader. Input files are never modified.

| Relative path | Required structure |
| --- | --- |
| `metadata/bkv_sample_metadata.csv` | `sample`, `group`; one row per source profile. Group labels must be exactly `control` or `treatment`. The historical `treatment` label denotes BKPyV-DNAemia case status, not random allocation. |
| `metadata/bkv_sample_metadata_16_qc.csv` | `sample`, `group`; historical retained subset and column order. The computed QC sample set must match it. |
| `mirdeep/counts_matrix.csv` | First column `miRNA`, unique human-feature IDs; remaining columns are sample IDs. Nonnegative finite human-miRNA counts, including features subsequently filtered. |
| `mirdeep/bkv_counts_matrix.csv` | First column `miRNA`; sample columns and two viral rows, `bkv-miR-B1-3p` and `bkv-miR-B1-5p`. Nonnegative finite counts. The annotated 3p sequence is shared with JCPyV. |
| `mirdeep/volcano_results_hybrid_human_qn_bkv_targeted_16_qc.csv` | Historical comparison export. Consumed columns: `miRNA`, `log2FC`, `pvalue`. Extra columns may remain. Used to assert reproduction of original transformed mean differences and p-values. |
| `metaphlan_output/merged_abundance_table_species_phyloseq.csv` | First column `TAXA`; sample columns containing nonnegative finite bacterial relative percentages. All source profiles should sum to 100 within rounding tolerance (0.001 percentage points). |
| `metaphlan_output/merged_abundance_table_taxa_phyloseq.csv` | First column `TAXA`; matching row IDs and columns `Kingdom`, `Class`, `Family`, `Genus`, `Species`. Additional `Phylum` and `Order` columns are permitted. All `Kingdom` values must be `k__Bacteria`. Rank prefixes such as `c__` are stripped for display. |
| `mirdeep/correlation_miRNA_microbiome_all_classes_16_qc_9v7_rho.csv` | First column containing miRNA IDs; bacterial class names as columns, containing historical Spearman rho values. Used for numerical reproduction checks. |

All sample IDs must agree exactly across the tables; retain the historical IDs for the current script's numeric sample ordering in the Sankey. Taxonomic IDs must match between abundance and taxonomy tables. Missing taxonomy ranks are grouped as unclassified labels; missing abundance/count values are invalid.

The original metadata had 22 profiles (11 cases and 11 controls). Human-miRNA QC retained 16 (9 cases and 7 controls). These counts do not establish correspondence with the earlier published 8-versus-8 cohort. A clinical sample crosswalk remains necessary.

The workflow starts from quantified profiles and cannot recover upstream sequencing QC, residual human DNA, batches, negative controls, viruses/fungi/parasites absent from these tables, or calibrated concentrations. Do not substitute normalized read yields for absolute microbial burden.

Every execution records input SHA-256 hashes and analysis package versions in its local `results/run_summary.json`. This manifest is intentionally not committed because local outputs are outside the code repository.
