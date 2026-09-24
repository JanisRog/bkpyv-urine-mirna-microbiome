> This GitHub directory contains code only. The local author package described below additionally includes private inputs, figures, results and manuscript drafts. Obtain approved inputs separately. No original sequencing inputs or manuscript text are tracked here.

# BKPyV reanalysis for author review

This is the current analysis code as of 24 September 2026. It supersedes the earlier quantile-normalized analysis and the older merged bacterial table for the revised results. Previous outputs remain preserved outside this package. Author review and source-metadata resolution remain necessary before journal submission.

The revised human counts remove repeated assignments across precursors and exclude ambiguity across distinct mature names. Main analysis: at least 2,000 assigned human counts, 9 cases/7 controls. Sensitivities: all 22, at least 1,000 counts, previous combined criterion, and individual influence exclusions. DESeq2 and edgeR test integer counts. Bacterial profiles were rebuilt from 22 accession-level MetaPhlAn outputs with verified archive mapping, yielding 119 species labels. The all-22 bacterial analysis is primary for that component. See `ANALYSIS_PLAN.md` for retrospective decisions and limitations.

No human miRNA survives FDR in any of the four specified inclusion strategies with either count model. The S24-excluded influence analysis does yield discoveries and must not be mistaken for a justified exclusion. Bacterial and group-conditioned integrated tests do not survive their stated FDR corrections. Viral results depend on normalization. Absence of corrected findings does not establish equivalence.

The 22-sample Kraken2 screen and competitive BKPyV/JCPyV DNA alignment are complete. All 22 urine DNA libraries, including blood-DNAemia-negative controls, had deep BKPyV genome coverage. S14 was the only sample with deep JCPyV genome coverage. Other low-level JCPyV signals cannot be called definitively absent. The shared BK/JCPyV 3p miRNA remains non-specific to either virus. Sequence reads do not establish urine copies/mL or productive infection.

## Contents

- Manuscript and response drafts remain in the private local author package, not this code repository.
- `figures/`: synchronized figures in the private local author package, including the retained descriptive Sankey.
- `results/count_models/`, `results/tables/`, `results/dna_qc/`: model results, full taxonomic descriptions, sensitivities and DNA provenance checks.
- `inputs/`: local copies of small source tables and QC files; no FASTQ or BAM files. These inputs are not included in GitHub.
- `verified_data/`: bacterial tables rebuilt from accession profiles and the accession-to-sample mapping.
- `analysis/`: reproducible numerical analysis; the manuscript builder additionally depends on the prior revision's formatting and source text and is not part of the standalone analysis command.
- `cluster/`: count reconstruction, DNA inventory, 22-sample Kraken2 screening and completed competitive BKPyV/JCPyV alignment jobs. `run_dna_kraken2.sbatch` uses the existing PlusPF reference and host-filtered FASTQs. `run_polyoma_alignment.sbatch` uses the two fixed NCBI RefSeq genomes and reports high-confidence alignment counts and breadth.
- `analysis/summarize_kraken2.py`: validates the compact archives from all 22 array tasks against their run metadata and known sample mapping, then creates per-sample domain counts and the complete nonzero taxon table. It does not treat classifier assignments as verified infections or absolute abundance.
- `analysis/summarize_polyoma_reports.py`: validates all 22 per-sample BKPyV/JCPyV alignment archives, recomputes coverage breadth from their depth files, and joins the results with the Kraken2 screen.
- `analysis/final_dna_integration.py`: exports sample-level viral DNA, viral-miRNA and QC values and makes the polyomavirus and bacterial individual-level figures. This integration is exploratory and follows the prespecified table analysis.
- `analysis/compile_eukaryote_screen.py`: exports fungal and other non-metazoan eukaryotic species labels as a descriptive screen, without infection calls.
- `analysis/build_final_figure_pack.py`: assembles Figures 1–4 and retains the descriptive Sankey as Figure 4, with its interactive HTML version.

## Rerun the numerical analysis

Use Python 3.12 with `pip install -r analysis/requirements.txt`. The recorded run used R 4.5.1, DESeq2 1.50.2 and edgeR 4.8.2; full R session information is in `results/count_models/R_sessionInfo.txt`. Install compatible Bioconductor packages in an isolated library before running. No dependency installation occurs automatically. The local `Rlib` folder is not included in the archive or repository.

From this package directory, with the appropriate Python environment active:

```bash
python analysis/run_revision.py \
  --profiles inputs/dna_profiles \
  --mapping results/dna_qc/ena_sample_mapping.tsv \
  --metadata verified_data/metadata/bkv_sample_metadata.csv \
  --recount inputs/recount \
  --qc inputs/smallRNA_sample_qc.csv \
  --rlib Rlib \
  --out rerun_output
```

The output folder must not already exist. If packages are installed in a standard R library, supply an existing empty directory to `--rlib`; standard R libraries remain on the search path. Use `--rscript` if Rscript is not on PATH. On this Mac, the existing package-local Rlib can be used directly. The notebook `analysis/reanalysis_v2.ipynb` provides the same run with explicit path prompts.

`run_revision.py` rebuilds bacterial tables, fits the six count-model configurations, runs permutation tests, creates figures and independently checks saved results. The completed 24 September rerun reproduced the earlier key tables byte for byte. It does not redo raw-read alignment or generate manuscript text. The direct validator is:

```bash
python analysis/validate_results.py --root . --recount inputs/recount
python -m unittest discover -s tests -v
```

## Work remaining before submission

1. Verify original tool/reference releases, DNA trimming history, library metrics and batches where source records exist. Current container versions alone will not establish historical versions.
2. Resolve urine fraction, input volume, the E. coli UTI control's identity, individual viral PCR/pathology data, collection timing and the original six miRNA exclusions with the source investigators. If unavailable, retain explicit limitations. No absolute abundance per mL can be recovered from the current relative tables alone.
3. Author review of conclusions, figure choices, article-length limits, citations and all response statuses; then an approved public code/archive release and data-access statement. The repository remains private.

Full source input and output hashes are recorded in `package_checksums.json`. Clinical sequence inputs and manuscript drafts must not be committed to GitHub. A private local archive containing the small input tables is for authorized collaborators only.
