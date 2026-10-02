# BKPyV urinary miRNA and microbiome analysis

This page preserves the earlier quantile-normalized workflow for audit. It does not generate the results or figure numbering in the revised brief communication; see [`reanalysis_v2/README.md`](../reanalysis_v2/README.md) for the current workflow.

The workflow reproduces the submitted table-based analysis, reports normalization sensitivities, summarizes bacteria from class to species, and produces the descriptive Sankey and supplementary correlation heatmaps. It does not process raw sequencing reads or establish validated biomarkers, biological interactions or absolute microbial loads.

## Install

Use Python 3.12. From the repository folder, create an isolated environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r analysis/requirements.txt
```

On Windows, activate with `.venv\Scripts\activate` instead. Analysis dependencies are pinned to the versions used in the recorded revision run.

## Run from the command line

Point `--data-dir` at the existing input folder described in [the historical input list](INPUTS.md):

```bash
python analysis/reanalyze.py --data-dir /path/to/approved/input_data --out-dir outputs
```

Replace the example input path. Outputs must be outside the input folder. The command writes to `outputs/results/` and `outputs/figures/`; rerunning with the same output folder replaces those generated files. Inputs are read only.

The original inputs and historical exported comparison tables are required. This is a reproduction workflow for this study, not a drop-in analysis for an arbitrary new cohort: QC selection and historical effects/correlations are checked against the supplied exports. Those checks intentionally fail if incompatible inputs are used.

## Run in VS Code or Jupyter

1. Open the repository folder in VS Code, install the Python and Jupyter extensions if needed, and select `.venv` as the notebook kernel.
2. With the environment active, install notebook support: `python -m pip install jupyterlab ipykernel`.
3. Open `analysis/revision_analysis.ipynb` and run its cells. Enter your input folder when prompted. No personal path needs to be saved in the notebook source.

Alternatively, launch `python -m jupyterlab analysis/revision_analysis.ipynb` from the repository folder. You can set the `BKPYV_DATA_DIR` environment variable before launching instead of using the input prompt. The notebook finds the repository from its working directory and saves outputs under `outputs/`.

## Interpretation

- Human miRNAs use log2 transformation followed by the submitted minimum-rank quantile normalization. Viral miRNAs use unnormalized log2 counts in the reproduced analysis.
- Averaged tied quantiles and viral counts per million human-miRNA counts are separate sensitivities. The human-miRNA denominator is not total sequencing depth or an absolute viral-load estimate.
- The submitted candidate set is retained for traceability. The minimum-rank tie treatment is not endorsed as the final preferred model. Authors must justify the final primary normalization/model before submission.
- The revision reproduced 268 original miRNA tests and 168 correlations; none passed their stated global FDR correction. Class-to-species tests are exploratory percentage comparisons with within-rank and across-rank FDR.
- The Sankey connects measurements through sample identity. Its left and right link sets have independent units and scales. A display cutoff is not a significance threshold, and links do not establish interactions.

See [the historical analysis notes](ANALYSIS.md) for the statistical scope and output definitions.

## Data access and repository contents

No clinical records, sequencing inputs, manuscript drafts or generated sample-level results are included. Access to the miRNA inputs is subject to the original data-sharing permissions. Request approved inputs from the study authors; no public access route is asserted here. Code alone cannot reproduce the clinical results without those inputs.

The ignore rules exclude local data, outputs, environments and credentials. Before pushing notebook changes, clear all cell outputs and check that no personal path or data was saved. The automated repository check helps detect common mistakes but does not replace inspection of `git diff --cached`.

## Checks and publication

```bash
python -m unittest discover -s tests -v
python tests/check_repository.py
```

CI runs these checks on synthetic numerical examples without accessing study data. It does not independently validate the clinical results. Local full-data reproduction remains a separate check.

No reuse license or archived release DOI has been assigned yet. See [the release checklist](RELEASE_CHECKLIST.md) before citing a frozen release. Code assistance and author responsibility should be disclosed consistently with the manuscript.
