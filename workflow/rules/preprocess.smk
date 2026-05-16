rule preprocess_impute:
    input:
        feature_path=feature_input,
        metadata_path=metadata_input,
    output:
        processed=f"{RESULTS_DIR}/preprocess_impute/{{dataset_id}}/bulk_x_features.csv",
        benchmark=f"{RESULTS_DIR}/preprocess_impute/{{dataset_id}}/benchmark.rds",
        summary=f"{RESULTS_DIR}/preprocess_impute/{{dataset_id}}/preprocess_summary.csv",
    log:
        f"{LOGS_DIR}/preprocess_impute/{{dataset_id}}.log",
    params:
        missingness_threshold=config["preprocess"]["missingness_threshold"],
        min_unique_values=config["preprocess"]["min_unique_values"],
        impute=bool_str(config["preprocess"]["impute"]),
        exclude_missing_samples=bool_str(config["preprocess"].get("exclude_missing_samples", False)),
        run_imputation_benchmark=bool_str(config["preprocess"]["run_imputation_benchmark"]),
        use_benchmark_for_feature_filtering=bool_str(
            config["preprocess"].get("use_benchmark_for_feature_filtering", True)
        ),
        run_distribution_shift_filter=bool_str(config["preprocess"]["run_distribution_shift_filter"]),
        benchmark_methods=join_comma(config["preprocess"]["benchmark_methods"]),
        imputation_method=config["preprocess"]["imputation_method"],
        imputation_method_arg=(
            ""
            if config["preprocess"]["imputation_method"] is None
            else f"--imputation_method {config['preprocess']['imputation_method']}"
        ),
        na_frequencies=join_comma(config["preprocess"]["na_frequencies"]),
        error_threshold=config["preprocess"]["error_threshold"],
        error_metric=config["preprocess"]["error_metric"],
        distribution_shift_z_threshold=config["preprocess"]["distribution_shift_z_threshold"],
        distribution_shift_ks_alpha=config["preprocess"].get("distribution_shift_ks_alpha", 0.05),
        distribution_shift_group_key_arg=(
            ""
            if config["preprocess"].get("distribution_shift_group_key") is None
            else f"--distribution_shift_group_key {config['preprocess']['distribution_shift_group_key']}"
        ),
        min_unique_values_arg=(
            ""
            if config["preprocess"]["min_unique_values"] is None
            else f"--min_unique_values {config['preprocess']['min_unique_values']}"
        ),
    shell:
        r"""
        mkdir -p "$(dirname {output.processed:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/preprocess_impute.R \
          --feature_path {input.feature_path:q} \
          --metadata_path {input.metadata_path:q} \
          --sample_key {SAMPLE_KEY:q} \
          --out_dir "$(dirname {output.processed:q})" \
          --missingness_threshold {params.missingness_threshold} \
          {params.min_unique_values_arg} \
          --impute {params.impute:q} \
          --exclude_missing_samples {params.exclude_missing_samples:q} \
          --run_imputation_benchmark {params.run_imputation_benchmark:q} \
          --use_benchmark_for_feature_filtering {params.use_benchmark_for_feature_filtering:q} \
          --benchmark_methods {params.benchmark_methods:q} \
          {params.imputation_method_arg} \
          --na_frequencies {params.na_frequencies:q} \
          --error_threshold {params.error_threshold} \
          --error_metric {params.error_metric:q} \
          --run_distribution_shift_filter {params.run_distribution_shift_filter:q} \
          --distribution_shift_z_threshold {params.distribution_shift_z_threshold} \
          --distribution_shift_ks_alpha {params.distribution_shift_ks_alpha} \
          {params.distribution_shift_group_key_arg} \
          > {log:q} 2>&1
        """


rule evaluate_preprocessing:
    input:
        feature_paths=lambda wildcards: [feature_input_path(dataset_id) for dataset_id in DATASET_IDS],
        benchmark_paths=expand(
            f"{RESULTS_DIR}/preprocess_impute/{{dataset_id}}/benchmark.rds",
            dataset_id=DATASET_IDS
        ),
        preprocess_summary_paths=expand(
            f"{RESULTS_DIR}/preprocess_impute/{{dataset_id}}/preprocess_summary.csv",
            dataset_id=DATASET_IDS
        ),
    output:
        plot=f"{EVALUATE_PREPROCESSING_ROOT}/error_distribution.pdf",
        error_threshold=f"{EVALUATE_PREPROCESSING_ROOT}/error_threshold.pdf",
        error_threshold_detailed=f"{EVALUATE_PREPROCESSING_ROOT}/error_threshold_detailed.pdf",
        error_threshold_log10n=f"{EVALUATE_PREPROCESSING_ROOT}/error_threshold_log10n.pdf",
        error_threshold_log10n_detailed=f"{EVALUATE_PREPROCESSING_ROOT}/error_threshold_log10n_detailed.pdf",
        retention_summary_csv=f"{EVALUATE_PREPROCESSING_ROOT}/feature_retention_summary.csv",
        retention_summary_pdf=f"{EVALUATE_PREPROCESSING_ROOT}/feature_retention_summary.pdf",
        retention_summary_detailed_pdf=f"{EVALUATE_PREPROCESSING_ROOT}/feature_retention_summary_detailed.pdf",
        retention_summary_log10n_pdf=f"{EVALUATE_PREPROCESSING_ROOT}/feature_retention_summary_log10n.pdf",
        retention_summary_log10n_detailed_pdf=f"{EVALUATE_PREPROCESSING_ROOT}/feature_retention_summary_log10n_detailed.pdf",
    threads: RULE_THREADS.get("evaluate_preprocessing", 1)
    resources:
        mem_mb=rule_resource("evaluate_preprocessing", "mem_mb", 16000),
    log:
        f"{LOGS_DIR}/evaluate_preprocessing.aggregate.{config['preprocess']['error_metric']}.log",
    params:
        error_metric=config["preprocess"]["error_metric"],
        imputation_method=config["preprocess"]["imputation_method"],
        imputation_method_arg=(
            ""
            if config["preprocess"]["imputation_method"] is None
            else f"--imputation_method {config['preprocess']['imputation_method']}"
        ),
        thresholds=join_comma(config["plot_benchmark"]["thresholds"]),
        feature_table_dir=EVALUATE_FEATURE_TABLE_DIR,
        feature_filename=EVALUATE_FEATURE_FILENAME,
        benchmark_dir=f"{RESULTS_DIR}/preprocess_impute",
    shell:
        r"""
        mkdir -p "$(dirname {output.plot:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/evaluate_preprocessing.R \
          --feature_table_dir {params.feature_table_dir:q} \
          --feature_filename {params.feature_filename:q} \
          --benchmark_dir {params.benchmark_dir:q} \
          --error_metric {params.error_metric:q} \
          {params.imputation_method_arg} \
          --n_jobs {threads} \
          --thresholds {params.thresholds:q} \
          > {log:q} 2>&1
        """
