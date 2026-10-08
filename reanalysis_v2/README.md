> This GitHub directory contains code only. The local author package described below additionally includes private inputs, figures, results and manuscript drafts. Obtain approved inputs separately. No original sequencing inputs or manuscript text are tracked here.

# BKPyV reanalysis for author review

This is the analysis code as of 8 October 2026. It supersedes the earlier quantile-normalized analysis and the older merged bacterial table for the revised results. The study independently analyses previously collected sequencing data to ask whether urinary bacterial composition is associated with BKPyV reactivation. Previous outputs remain preserved outside this package. Author review and source-metadata resolution remain necessary before journal submission.

The revised human counts remove repeated assignments across precursors and exclude ambiguity across distinct mature names. Main analysis: at least 2,000 assigned human counts, 9 cases/7 controls. Sensitivities: all 22, at least 1,000 counts, previous combined criterion, and individual influence exclusions. DESeq2 and edgeR test integer counts. Bacterial profiles were rebuilt from 22 accession-level MetaPhlAn outputs with verified archive mapping, yielding 119 species labels. The all-22 bacterial analysis is primary for that component. See `ANALYSIS_PLAN.md` for retrospective decisions and limitations.

Main Figure 1A places human and viral miRNAs in one DESeq2 negative-binomial model, with size factors estimated from the human counts before adding the two viral rows. Fifteen human miRNAs meet the nominal display rule in the main 9/7 subset (11 higher and four lower in cases); none has q < .05. The shared BKPyV/JCPyV 3p sequence has q=0.057, while BKPyV 5p has no primary p-value because S7 triggers Cook's-distance filtering. In the 11/11 sensitivity, BKPyV 5p has DESeq2 q=0.012 but robust edgeR q=0.374, so the result is method- and subset-sensitive. The S24-excluded human influence analysis yields a discovery but does not justify excluding S24. Figure 1B shows all 22 MetaPhlAn bacterial class profiles and equal-weight group means. Main Figure 2 is a descriptive Sankey that links class and miRNA patterns by sample identity; its ribbons do not test associations. The host-normalized, group-adjusted miRNA–genus heatmap is Supplementary Figure S3. None of its 221 pairs survives FDR correction (minimum q=0.876). Bacterial and integrated results remain exploratory.

The 22-sample Kraken2 screen and competitive BKPyV/JCPyV DNA alignment are complete. All 22 urine DNA libraries, including blood-DNAemia-negative controls, had deep BKPyV genome coverage. S14 was the only sample with deep JCPyV genome coverage. Other low-level JCPyV signals cannot be called definitively absent. The shared BK/JCPyV 3p miRNA remains non-specific to either virus. Sequence reads do not establish urine copies/mL or productive infection.

## Contents

- Manuscript and response drafts remain in the private local author package, not this code repository.
- `figures/`: synchronized figures in the private local author package, including the main Sankey and supplementary miRNA–genus heatmap.
- `results/count_models/`, `results/tables/`, `results/dna_qc/`: model results, full taxonomic descriptions, sensitivities and DNA provenance checks.
- `inputs/`: local copies of small source tables and QC files; no FASTQ or BAM files. These inputs are not included in GitHub.
- `verified_data/`: bacterial tables rebuilt from accession profiles and the accession-to-sample mapping.
- `analysis/`: reproducible numerical analysis; the manuscript builder additionally depends on the prior revision's formatting and source text and is not part of the standalone analysis command.
- `cluster/`: count reconstruction, DNA inventory, 22-sample Kraken2 screening and completed competitive BKPyV/JCPyV alignment jobs. `run_dna_kraken2.sbatch` uses the existing PlusPF reference and host-filtered FASTQs. `run_polyoma_alignment.sbatch` uses the two fixed NCBI RefSeq genomes and reports high-confidence alignment counts and breadth.
- `analysis/summarize_kraken2.py`: validates the compact archives from all 22 array tasks against their run metadata and known sample mapping, then creates per-sample domain counts and the complete nonzero taxon table. It does not treat classifier assignments as verified infections or absolute abundance.
- `analysis/summarize_polyoma_reports.py`: validates all 22 per-sample BKPyV/JCPyV alignment archives, recomputes coverage breadth from their depth files, and joins the results with the Kraken2 screen.
- `analysis/final_dna_integration.py`: exports sample-level viral DNA, viral-miRNA and QC values and makes the polyomavirus and bacterial individual-level figures. This integration is exploratory and follows the prespecified table analysis.
- `analysis/compile_eukaryote_screen.py`: exports fungal and other non-metazoan eukaryotic species labels as a descriptive screen, without infection calls.
- `analysis/run_joint_mirna_models.R`: fits 244 human features and two prespecified viral sequences together across six subsets, with fixed human-derived normalization in DESeq2 and robust edgeR. It saves Cook's-distance and S7 influence diagnostics. The BK/JC shared 3p sequence remains species-ambiguous.
- `analysis/make_joint_mirna_figure.py`: makes the earlier single-panel miRNA figure from the primary joint DESeq2 fit. BKPyV 5p is shown in a separate untested band rather than being assigned an invented volcano p-value.
- `analysis/make_owner_focus_figures.py`: makes the current Figure 1A/B and main Figure 2 from the verified private model, recount, metadata and bacterial tables. It writes the figure captions and Sankey link table. Input data are excluded from GitHub.
- `analysis/make_final_human_figure.py`: preserves the previous human-only Figure 1 generation for audit; it is superseded by the joint figure above.
- `analysis/make_mirna_genus_heatmap.py`: makes the single-panel heatmap now used as Supplementary Figure S3 from host-normalized miRNA ranks and genus centred-log-ratio ranks. It removes case/control means, uses 19,999 within-group permutations, corrects all 221 pairs together, and exports raw-count, all-22 and pseudocount sensitivities. Candidate selection and correlation use the same cohort, so these are exploratory tests.
- `analysis/build_tid_brief_figures.py`: preserves the 2 October layout, with the heatmap as main Figure 2 and Sankey as Supplementary Figure S3. It uses explicit input directories and checks the 22-sample viral table; it does not assemble the 8 October figure arrangement.
- `analysis/build_three_figure_pack.py` and `analysis/build_final_figure_pack.py`: preserve earlier layouts for audit; they do not represent the final brief communication.

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

`run_revision.py` rebuilds bacterial tables, fits the six count-model configurations, runs permutation tests, creates figures and independently checks saved results. The completed 24 September rerun reproduced the earlier key tables byte for byte. It does not redo raw-read alignment or generate manuscript text.

To generate the earlier human-only figure from a completed run, use a fresh output directory:

```bash
python analysis/make_final_human_figure.py \
  --models rerun_output/results/count_models \
  --out human_figure_original_style
```

After `run_revision.py` has generated host-only normalization factors, keep the private recount input in the ignored rerun directory and run the joint models and association heatmap there:

```bash
mkdir -p rerun_output/inputs
cp -R inputs/recount rerun_output/inputs/recount
Rscript analysis/run_joint_mirna_models.R rerun_output rerun_output/joint_mirna_models Rlib
python analysis/make_joint_mirna_figure.py \
  --models rerun_output/joint_mirna_models \
  --out rerun_output/joint_figure_1
python analysis/make_mirna_genus_heatmap.py \
  --root rerun_output --models rerun_output/joint_mirna_models \
  --out rerun_output/mirna_genus_heatmap
```

The output directories must not already exist. The joint model reads the local recount matrices and the host-only factors saved by the preceding run. The following command assembles the archived 2 October layout, with the association heatmap as Figure 2. Supply the independently verified DNA-integration output:

```bash
python analysis/build_tid_brief_figures.py \
  --root rerun_output \
  --joint-figure rerun_output/joint_figure_1 \
  --heatmap rerun_output/mirna_genus_heatmap \
  --integration /path/to/validated_dna_integration \
  --out rerun_output/tid_figures
```

The 8 October Figure 1A/B and Figure 2 are generated by `analysis/make_owner_focus_figures.py` from the same verified private outputs. The script expects `joint_mirna_models_2026-10-02_v2/main_ge2000`, `verified_data/metadata/bkv_sample_metadata.csv`, `results/tables/class_profiles_all22.csv`, and `inputs/recount` under the local `reanalysis_v2` package root. From that package, run `python analysis/make_owner_focus_figures.py`; it writes to `project_owner_focus_2026-10-08/figures`. These files are excluded from GitHub. The DNA-integration input is generated by `analysis/final_dna_integration.py` from the verified Kraken2 and competitive BKPyV/JCPyV alignment summaries; those large upstream reads are not in this repository. To check a numerical run:

```bash
python analysis/validate_results.py --root . --recount inputs/recount
python -m unittest discover -s tests -v
```

## Work remaining before submission

1. Verify original tool/reference releases, DNA trimming history, library metrics and batches where source records exist. Current container versions alone will not establish historical versions.
2. Resolve urine fraction, input volume, the E. coli UTI control's identity, individual viral PCR/pathology data, collection timing and the original six miRNA exclusions with the source investigators. If unavailable, retain explicit limitations. No absolute abundance per mL can be recovered from the current relative tables alone.
3. Author review of conclusions, article-length limits, citations and the response letter. Update the manuscript's code-availability statement with the public commit or release URL. A reuse license and archival DOI are separate decisions; neither is assigned here.

Full source input and output hashes are recorded in `package_checksums.json`. Clinical sequence inputs and manuscript drafts must not be committed to GitHub. A private local archive containing the small input tables is for authorized collaborators only.
