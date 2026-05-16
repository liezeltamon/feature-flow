# Snakemake workflow: `feature-flow`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/liezeltamon/feature-flow/actions/workflows/main.yaml/badge.svg?branch=main)](https://github.com/liezeltamon/feature-flow/actions/workflows/main.yaml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/liezeltamon/feature-flow)

A Snakemake workflow for reducing high-dimensional feature tables into module-level summaries with statistical testing and optional biological annotation.

- [Snakemake workflow: `feature-flow`](#snakemake-workflow-feature-flow)
  - [Overview](#overview)
  - [Input data](#input-data)
  - [Output](#output)
  - [Usage](#usage)
  - [Deployment options](#deployment-options)
  - [Workflow profiles](#workflow-profiles)
  - [Authors](#authors)
  - [References](#references)

## Overview

The workflow is built using [Snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. **Input combination** *(optional)* — Merges multiple feature tables via row-bind (same features, different samples) or column-bind (same samples, different features) before downstream processing.
2. **Feature normality testing** — Tests each feature for normality (Shapiro-Wilk) and records results as an annotation source for downstream summary selection.
3. **Feature preprocessing and imputation** — Filters features by missingness and unique-value thresholds, optionally imputes missing values, optionally benchmarks imputation methods, and optionally applies distribution-shift filtering.
4. **Module reduction** — Reduces preprocessed features into hierarchical-clustering modules using a configured discovery subset, with optional repeated subsampling when imputation is applied.
5. **Module annotation** *(optional)* — Annotates gene-like modules with GO:BP pathway labels via `gprofiler2`, or names modules from top-loading members using an LLM.
6. **Module set enrichment** *(optional)* — Tests module members for enrichment in configured GMT gene sets.
7. **Univariate testing — modules** — Fits per-module models (e.g. `lmer`, rank-transform) and stores effect sizes, p-values, and model-assumption results.
8. **Univariate testing — original features** — Same testing on original/member features, dropping NA samples per feature and requiring a minimum group size.
9. **Result summarisation** — Converts method-specific result matrices into wide summary CSVs, selecting parametric or non-parametric fields based on model assumptions met.
10. **Visualisation** — Plots module-level summaries and original feature distributions with group labels and optional module annotations.

The workflow supports three input modes: **single** (one feature table), **batch** (multiple dataset directories discovered automatically), and **combine** (tables merged before processing).

Detailed information about input data and workflow configuration can be found in the [`config/README.md`](config/README.md).

## Input data

| Input | Description |
| --- | --- |
| Feature table | CSV with `sample_key` as the first column; all subsequent numeric columns are treated as features |
| Sample metadata | CSV with `sample_key` column plus grouping, individual ID, batch, and covariate columns |
| Dataset directories *(batch mode)* | Subdirectories under `input_root`, each containing `bulk_x_features.csv` and `bulk_metadata.csv` |
| GMT file *(optional)* | Gene set file for module set enrichment |

## Output

All outputs are written to `<results_dir>/` as configured in `config/config.yaml`.

| Directory | Key output files |
| --- | --- |
| `preprocess_impute/{dataset_id}/` | `bulk_x_features.csv`, `benchmark.rds`, `preprocess_summary.csv` |
| `evaluate_preprocessing/{error_metric}/` | `error_distribution.pdf`, `error_threshold*.pdf`, `feature_retention_summary*.csv/pdf` |
| `modules_hc/{dataset_id}/{subset_name}/` | `bulk_x_features.csv`, `module_eigengenes_hc.rds`, `module_members/` |
| `meaning_module_genes/{dataset_id}/{subset_name}/` | `bulk_x_features.csv`, `module_annotations.csv`, `modules/` |
| `meaning_module_llm/{dataset_id}/{subset_name}/` | `bulk_x_features.csv`, `module_annotations.csv`, `modules/` |
| `univariate_test/{dataset_id}/{method}/` | `result_matrices.rds`, `failed_features.csv` |
| `summarise_tests/{dataset_id}/` | `summary_df.csv` |
| `summarise_tests/` | `summary_df.csv` (aggregate across all datasets) |
| `univariate_test_members/{dataset_id}/{method}/` | `result_matrices.rds`, `failed_features.csv` |
| `summarise_tests_members/{dataset_id}/` | `summary_df.csv` |
| `summarise_tests_members/` | `summary_df.csv` (aggregate across all datasets) |

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/liezeltamon/feature-flow).

If you use this workflow in a paper, please cite the repository URL or its DOI and the tools listed in the [References](#references) section.

## Deployment options

Change to the workflow directory and adjust options in `config/config.yaml`.

```bash
cd path/to/feature-flow
```

Perform a dry run to check the workflow before execution:

```bash
snakemake --dry-run
```

Run with test files using **conda**:

```bash
snakemake --cores 2 --sdm conda --directory .test
```

Run with **apptainer** / **singularity**:

```bash
snakemake --cores 2 --sdm conda apptainer --directory .test
```

Run on an HPC cluster via **SLURM** (recommended for production). First activate the environment:

```bash
bash envs/feature-module-reduction-env.sh
load_mamba
mamba activate feature-module-reduction-env
pip install snakemake-executor-plugin-slurm
```

Then submit using the provided template script:

```bash
bash slurm_submit_template.sh
```

For runs that use LLM module naming, source OpenAI credentials before submitting:

```bash
source ~/.config/openai/env.sh
```

## Workflow profiles

The `profiles/` directory can contain any number of [workflow-specific profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) that users can choose from.
The [profiles `README.md`](profiles/README.md) provides more details.

## Authors

- Liezel Tamon
  - University of Oxford
  - [ORCID profile](https://orcid.org/0000-0003-3705-6019)

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2

> Bashford-Rogers Lab. _vdjremix_. https://github.com/Bashford-Rogers-lab/vdjremix
