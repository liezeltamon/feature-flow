rule summarise_tests_members:
    input:
        results=expand(
            f"{RESULTS_DIR}/univariate_test_members/{{dataset_id}}/{{method}}/result_matrices.rds",
            dataset_id=DATASET_IDS,
            method=config["univariate"]["methods"]
        ),
        normality_results=expand(
            f"{RESULTS_DIR}/test_normality/{{dataset_id}}/normality_test_results.csv",
            dataset_id=DATASET_IDS
        ),
    output:
        summary=f"{RESULTS_DIR}/summarise_tests_members/summary_df.csv",
        per_source=expand(
            f"{RESULTS_DIR}/summarise_tests_members/{{dataset_id}}/summary_df.csv",
            dataset_id=DATASET_IDS
        ),
    log:
        f"{LOGS_DIR}/summarise_tests_members.log",
    params:
        src_dir=f"{RESULTS_DIR}/univariate_test_members",
        normality_dir=f"{RESULTS_DIR}/test_normality",
        parametric_method=config.get("summarise_tests", {}).get("parametric_method", "lmer"),
        nonparametric_method=config.get("summarise_tests", {}).get("nonparametric_method", "rank_transform"),
    shell:
        r"""
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/summarise_tests.R \
          --src_dir {params.src_dir:q} \
          --module_normality_dir {params.normality_dir:q} \
          --out_dir "$(dirname {output.summary:q})" \
          --parametric_method {params.parametric_method:q} \
          --nonparametric_method {params.nonparametric_method:q} \
          > {log:q} 2>&1
        """
