# BKPyV urinary miRNA and microbiome analysis

Reproducible exploratory secondary analysis in kidney transplant recipients. The current revision uses corrected mature-miRNA integer counts, DESeq2 and edgeR, bacterial profiles rebuilt from accession-level outputs, compositional sensitivity analyses and a descriptive Sankey.

Start with [the current workflow](reanalysis_v2/README.md) and [the retrospective analysis specification](reanalysis_v2/ANALYSIS_PLAN.md). Open `reanalysis_v2/analysis/reanalysis_v2.ipynb` for an interactive run, or use `run_revision.py` in the same directory. The notebook prompts for actual approved input locations; no personal paths or execution outputs are stored in its source.

The previous scripts under `analysis/` reproduce the historical quantile-normalized analysis and are retained for traceability. They do not generate the current revised results. The [legacy workflow](docs/LEGACY_WORKFLOW.md) describes them.

## Scope and limitations

The current work is an author-review revision, not a completed or validated clinical pipeline. Further DNA viral/fungal/microbial-eukaryotic classification, residual-host assessment and clinical/laboratory provenance checks remain outstanding. Sample thresholds are retrospective and unvalidated. No subset is selected for producing significance; influence analyses are reported separately. Relative profiles and normalized miRNA counts are not absolute microbial or viral loads.

This private repository contains code and documentation only. No study input tables, clinical records, sequencing reads, manuscript drafts or generated sample-level results are tracked. Reproduction requires approved inputs from the authors. The local package additionally contains small input tables, numerical results and drafts; it is not the GitHub distribution.

## Checks

```bash
python -m pip install -r reanalysis_v2/analysis/requirements.txt
python -m unittest discover -s tests -v
python -m unittest discover -s reanalysis_v2/tests -v
python tests/check_repository.py
```

CI uses synthetic examples and source checks. Full-data model results are separately checked by `reanalysis_v2/analysis/validate_results.py`; CI has no access to study data and does not independently fit the clinical models. R 4.5.1 with DESeq2 1.50.2 and edgeR 4.8.2 was used for the recorded local reanalysis.

Public release, reuse licensing and archival DOI remain subject to author approval. See [the release checklist](docs/RELEASE_CHECKLIST.md). Code assistance and author responsibility should be disclosed consistently with the manuscript.
