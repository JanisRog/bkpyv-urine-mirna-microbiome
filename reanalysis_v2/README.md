> This GitHub directory contains code only. The local author package described below additionally includes private inputs, figures, results and manuscript drafts. Obtain approved inputs separately. No original sequencing inputs or manuscript text are tracked here.

# BKPyV reanalysis for author review

This is the current analysis package as of 23 September 2026. It supersedes the earlier quantile-normalized analysis and the older merged bacterial table for the revised results. Previous outputs remain preserved outside this package. This package is not ready for journal submission.

The revised human counts remove repeated assignments across precursors and exclude ambiguity across distinct mature names. Main analysis: at least 2,000 assigned human counts, 9 cases/7 controls. Sensitivities: all 22, at least 1,000 counts, previous combined criterion, and individual influence exclusions. DESeq2 and edgeR test integer counts. Bacterial profiles were rebuilt from 22 accession-level MetaPhlAn outputs with verified archive mapping, yielding 119 species labels. The all-22 bacterial analysis is primary for that component. See `ANALYSIS_PLAN.md` for retrospective decisions and limitations.

No human miRNA survives FDR in any of the four specified inclusion strategies with either count model. The S24-excluded influence analysis does yield discoveries and must not be mistaken for a justified exclusion. Bacterial and group-conditioned integrated tests do not survive their stated FDR corrections. Viral results depend on normalization. Absence of corrected findings does not establish equivalence.

## Contents

- `BKPyV_reanalysis_manuscript.docx` and `BKPyV_reanalysis_response.docx`: updated drafts with pending issues stated.
- `figures/`: synchronized figures, including the retained descriptive Sankey.
- `results/count_models/`, `results/tables/`, `results/dna_qc/`: model results, full taxonomic descriptions, sensitivities and DNA provenance checks.
- `inputs/`: local copies of small source tables and QC files; no FASTQ or BAM files. These inputs are not included in GitHub.
- `verified_data/`: bacterial tables rebuilt from accession profiles and the accession-to-sample mapping.
- `analysis/`: reproducible numerical analysis; the manuscript builder additionally depends on the prior revision's formatting and source text and is not part of the standalone analysis command.
- `cluster/`: count reconstruction, the completed DNA inventory job, and a Kraken2 DNA pilot/array job. `run_dna_kraken2.sbatch` uses the existing PlusPF reference and host-filtered FASTQs. The first pilot completed but memory mapping took 16.4 hours for 1.39 million pairs. The revised script requests 128 GB to hold the database in RAM; rerun array index 3 and review its compact report before submitting all 22 samples.

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

`run_revision.py` rebuilds bacterial tables, fits the six count-model configurations, runs permutation tests, creates figures and independently checks saved results. It does not redo raw-read alignment, generate manuscript text, or complete the outstanding DNA classification. The direct validator is:

```bash
python analysis/validate_results.py --root . --recount inputs/recount
python -m unittest discover -s tests -v
```

## Work remaining before submission

1. Review the DNA Kraken2 pilot, then classify all 22 samples with the validated local PlusPF database. Audit viral, fungal and protozoan read support and residual-host classification; confirm BK/JC discrimination with nucleotide alignments/coverage before species-specific claims. The PlusPF convention does not guarantee exhaustive parasite coverage or identify this local database's exact build date.
2. Verify original tool/reference releases, DNA trimming history, library metrics and batches where records exist. Current container versions alone will not establish historical versions.
3. Resolve urine fraction, input volume, the E. coli UTI control's identity, individual viral PCR/pathology data, collection timing and the original six miRNA exclusions with the source investigators. Do not invent these facts; if unavailable, retain explicit limitations. No absolute abundance per mL can be recovered from the current relative tables alone.
4. Author review of conclusions, figure choices, article-length limits, citations and all response statuses; then an approved public code/archive release and data-access statement. The repository remains private.

Full source input and output hashes are recorded in `package_checksums.json`. Clinical sequence inputs and manuscript drafts must not be committed to GitHub. A private local archive containing the small input tables is for authorized collaborators only.
