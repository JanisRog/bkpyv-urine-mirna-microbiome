# BKPyV urinary miRNA and microbiome analysis

Reproducible exploratory analysis of previously collected urine sequencing data from kidney transplant recipients. The question is whether urinary bacterial composition is associated with BKPyV reactivation and co-varies with host or viral miRNAs. The current revision uses corrected mature-miRNA integer counts, DESeq2 and edgeR, bacterial profiles rebuilt from accession-level outputs, and compositional sensitivity analyses. The primary bacterial and miRNA comparisons use the same 16 recipients (nine cases, seven controls); the complete 22-sample data are retained for sensitivity analysis and in Supplementary Table S6. Main Figure 1 combines the miRNA volcano plot (A) and paired-cohort MetaPhlAn bacterial profiles (B); main Figure 2 is a descriptive Sankey. The miRNA–genus heatmap is Supplementary Figure S3.

Start with [the current workflow](reanalysis_v2/README.md) and [the retrospective analysis specification](reanalysis_v2/ANALYSIS_PLAN.md). Open `reanalysis_v2/analysis/reanalysis_v2.ipynb` for an interactive run, or use `run_revision.py` in the same directory. The notebook prompts for actual approved input locations; no personal paths or execution outputs are stored in its source.

The previous scripts under `analysis/` reproduce the historical quantile-normalized analysis and are retained for traceability. They do not generate the current revised results. The [legacy workflow](docs/LEGACY_WORKFLOW.md) describes them.

## Scope and limitations

The current work is an author-review revision, not a validated clinical pipeline. The 22-sample DNA viral and broad-domain screen is complete; source clinical/laboratory provenance checks remain outstanding. The human-miRNA figure uses the nominal descriptive display rule of unadjusted P < .05 and at least twofold change on corrected DESeq2 counts. This is not a rerun of the source CLC model; q values across all finite tests remain available. Sample thresholds are retrospective and unvalidated. No subset is selected for producing significance; influence analyses are reported separately. Relative profiles and normalized miRNA counts are not absolute microbial or viral loads.

This repository contains code and documentation only. No study input tables, clinical records, sequencing reads, manuscript drafts or generated sample-level results are tracked. Reproduction requires approved inputs from the authors. The local author package additionally contains small input tables, numerical results and drafts; it is not the GitHub distribution.

## Checks

```bash
python -m pip install -r reanalysis_v2/analysis/requirements.txt
python -m unittest discover -s tests -v
python -m unittest discover -s reanalysis_v2/tests -v
python tests/check_repository.py
```

CI uses synthetic examples and source checks. Full-data model results are separately checked by `reanalysis_v2/analysis/validate_results.py`; CI has no access to study data and does not independently fit the clinical models. R 4.5.1 with DESeq2 1.50.2 and edgeR 4.8.2 was used for the recorded local reanalysis.

No reuse license or archival DOI is assigned yet. Public visibility permits inspection of the code but does not change access to the study inputs. See [the release checklist](docs/RELEASE_CHECKLIST.md). Code assistance and author responsibility should be disclosed consistently with the manuscript.
