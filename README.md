#### Overview

- This directory contains a Snakemake pipeline for reducing large feature tables into module-level summaries and testing both modules and original member features.
- The pipeline supports:
  - single feature-table runs
  - batch runs over many feature-table directories
  - combined runs that merge multiple input tables first
- Main stages:
  - validate and optionally combine feature tables
  - test feature normality for annotation
  - preprocess, filter, and optionally impute features
  - reduce features into hierarchical-clustering modules
  - optionally annotate module meaning by pathway enrichment or LLM naming
  - run univariate tests on module features
  - run univariate tests on original/member features
  - export wide summary tables for downstream biology analyses

#### Environment

- R scripts run through the configured `rscript_bin` in the active config.
- Snakemake is required to orchestrate the workflow.
- Cluster execution uses the optional Snakemake SLURM executor plugin.
- LLM module naming requires OpenAI credentials when `module_meaning_llm.enabled: true`.
- Gene/pathway module naming uses `gprofiler2`, so compute nodes need network access to g:Profiler for that step.
- The repo-level `.Rprofile` and `renv::load()` are used by the R scripts; run from the repo root unless a script says otherwise.
- Minimal Snakemake environment setup:

```bash
bash envs/feature-module-reduction-env.sh
load_mamba
mamba activate feature-module-reduction-env
pip install snakemake-executor-plugin-slurm
```

#### Input data contract

- Feature tables must have the configured `sample_key` as the first column.
- Sample IDs stay as explicit columns, not row names.
- Numeric columns after `sample_key` are treated as features.
- Non-numeric feature-table columns after `sample_key` are ignored.
- Metadata tables provide grouping variables, individual IDs, batch variables, subset keys, and other model covariates.
- Feature-table sample IDs must be present in metadata; metadata is aligned to feature-table order.
- Batch mode expects immediate dataset directories containing:
  - `bulk_x_features.csv`
  - `bulk_metadata.csv`
- Combine mode supports:
  - `row_bind`: same features, different samples
  - `column_bind`: same samples, different features
- Contract-style feature IDs are kept intact as the official feature names.

#### Workflow scripts

- `Snakefile`
  - Defines the full workflow graph and final targets.
  - Supports `single`, `batch`, and `combine` input modes.
  - Discovers dataset IDs from input directories in batch mode.
  - Sends module-level testing to `univariate_test/` and `summarise_tests/`.
  - Sends original/member-feature testing to `univariate_test_members/` and `summarise_tests_members/`.
  - Chooses the canonical module table for downstream testing:
    - gene-like renamed module table when gene module meaning is enabled
    - LLM-renamed module table when LLM naming is enabled
    - raw `modules_hc` table otherwise

- `io_helpers.R`
  - Provides shared feature-table and metadata validation.
  - Preserves `sample_key` as character data so IDs are not silently converted.
  - Enforces sample-key presence and uniqueness.
  - Extracts numeric feature columns and records ignored non-numeric columns.
  - Aligns metadata rows to feature-table sample order.

- `combine_tables.R`
  - Combines multiple input datasets before downstream processing.
  - Uses `row_bind` when tables share feature columns and contain different samples.
  - Uses `column_bind` when tables share samples and contain different feature columns.
  - Writes combined feature table, combined metadata table, and an input manifest.

- `test_normality_wrapper.R`
  - Runs feature-level normality tests on numeric feature columns.
  - Writes `normality_test_results.csv` and `normality_summary.csv`.
  - Used as an annotation source for summaries.
  - Does not decide whether final univariate p-values come from parametric or nonparametric tests.

- `preprocess_impute.R`
  - Filters features by missingness and optional minimum unique values.
  - Optionally imputes missing values.
  - Optionally benchmarks imputation error across configured methods and missingness frequencies.
  - Optionally filters features by imputation benchmark performance.
  - Optionally applies distribution-shift filtering.
  - Writes processed `bulk_x_features.csv`, `benchmark.rds`, and `preprocess_summary.csv`.

- `evaluate_preprocessing.R`
  - Aggregates preprocessing benchmark outputs across datasets.
  - Produces imputation-error threshold plots.
  - Produces feature-retention summaries.
  - Uses the configured error metric, thresholds, and benchmark directory.

- `modules_hc.R`
  - Reduces preprocessed features into hierarchical-clustering modules.
  - Uses configured subset samples for module discovery.
  - Writes module eigengene tables, module member tables, and eigengene RDS output.
  - Uses repeated subsampling only when imputation was actually applied.
  - Can filter features by subset-specific unique values before module reduction.

- `meaning_module_genes.R`
  - Annotates gene-like modules using g:Profiler GO:BP enrichment.
  - Gene-like datasets are identified by dataset IDs containing `var_genes`.
  - Derives enrichment IDs from module member feature IDs using the configured separator.
  - Prefers significant pathway labels and falls back to suggestive labels when needed.
  - Writes renamed module tables, `module_annotations.csv`, and per-module audit files.

- `meaning_module_llm.py`
  - Names modules from top-loading module members using an LLM.
  - Intended for non-gene or mixed feature modules where pathway enrichment is not appropriate.
  - Writes renamed module table and `module_annotations.csv`.
  - Writes per-module prompt, raw response, parsed response, and meaning text for auditability.

- `module_set_enrichment.R`
  - Tests whether module members are enriched for configured GMT gene sets.
  - Supports module-member feature IDs with contract-style separators.
  - Can use only important-loading members or all members.
  - Writes `mhg_enrichment.csv` per dataset and subset.

- `univariate_test.R`
  - Fits per-feature models for each configured method, usually `lmer` and `rank_transform`.
  - Writes `result_matrices.rds` and `failed_features.csv`.
  - Stores effect-size, p-value, global p-value, assumption-test, and sample-count vectors.
  - Module mode uses `--missing_response_policy fail`.
  - Member mode uses `--missing_response_policy drop_feature_sample`.
  - Member mode drops NA samples separately per feature and requires at least 3 samples per requested group.
  - Feature values are scaled before testing, preserving existing pipeline behavior.
  - Uses an internal response column for modeling; output feature names remain unchanged.

- `summarise_tests.R`
  - Converts method-specific `result_matrices.rds` files into wide summary CSVs.
  - Writes one per-source `summary_df.csv` and one aggregate `summary_df.csv`.
  - Adjusts contrast p-values within each feature column.
  - Chooses final parametric vs nonparametric fields using `model_assumptions_met`.
  - `model_assumptions_met == TRUE` means final fields use the parametric method.
  - `model_assumptions_met == FALSE` means final fields use the nonparametric method.
  - `model_assumptions_met == NA` leaves final fields as `NA`.
  - Joins `is_normal` from normality outputs when available.

- `plot_features.R`
  - Plots module-level summaries and original feature distributions.
  - Uses metadata for group labels and sample covariates.
  - Uses module annotations when available.
  - Writes plot outputs under `plot_features/<dataset_id>/<subset_name>/`.

- `slurm_submit_template.sh`
  - Provides editable examples for launching Snakemake on SLURM.
  - Uses the configured run ID and config file.
  - Includes `--rerun-incomplete` and `--keep-going` for resumable cluster runs.
  - Sources OpenAI credentials for runs that use LLM module naming.

#### Output layout

- Preprocessing:

```text
<results_dir>/preprocess_impute/<dataset_id>/
  bulk_x_features.csv
  benchmark.rds
  preprocess_summary.csv
```

- Aggregate preprocessing evaluation:

```text
<results_dir>/evaluate_preprocessing/<error_metric>/
  error_distribution.pdf
  error_threshold*.pdf
  feature_retention_summary*.csv/pdf
```

- Module reduction:

```text
<results_dir>/modules_hc/<dataset_id>/<subset_name>/
  bulk_x_features.csv
  module_eigengenes_hc.rds
  module_members/
```

- Module meaning:

```text
<results_dir>/meaning_module_genes/<dataset_id>/<subset_name>/
  bulk_x_features.csv
  module_annotations.csv
  modules/

<results_dir>/meaning_module_llm/<dataset_id>/<subset_name>/
  bulk_x_features.csv
  module_annotations.csv
  modules/
```

- Module-level univariate testing:

```text
<results_dir>/univariate_test/<dataset_id>/<method>/
  result_matrices.rds
  failed_features.csv

<results_dir>/summarise_tests/<dataset_id>/
  summary_df.csv

<results_dir>/summarise_tests/
  summary_df.csv
```

- Member/original-feature univariate testing:

```text
<results_dir>/univariate_test_members/<dataset_id>/<method>/
  result_matrices.rds
  failed_features.csv

<results_dir>/summarise_tests_members/<dataset_id>/
  summary_df.csv

<results_dir>/summarise_tests_members/
  summary_df.csv
```

- Logs:

```text
<logs_dir>/<rule_name>/
<logs_dir>/<rule_name>.<details>.log
```

#### Running

- Run from the repo root.
- Use the default example config:

```bash
snakemake -s code/feature_module_reduction/Snakefile --cores 1
```

- Run the contract-style batch smoke test:

```bash
source ~/.config/openai/env.sh
snakemake -s code/feature_module_reduction/Snakefile \
  --configfile code/feature_module_reduction/config/examples/config_batch.yaml \
  --cores 1 -F
```

- Run the real `input_feature_table_level00` config:

```bash
snakemake -s code/feature_module_reduction/Snakefile \
  --configfile code/feature_module_reduction/config/input_feature_table_level00.yaml \
  --cores 1
```

- Run on SLURM using the template:

```bash
bash code/feature_module_reduction/slurm_submit_template.sh
```

- Main config templates:
  - `config/config.template.yaml`
  - `config/config_batch.template.yaml`
  - `config/config_combine_base.template.yaml`
  - `config/config_combine_row_bind.template.yaml`
  - `config/config_genes.template.yaml`
