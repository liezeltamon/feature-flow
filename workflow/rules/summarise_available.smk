rule summarise_tests_available:
    output:
        summary=f"{RESULTS_DIR}/summarise_tests_available/summary_df.csv",
        status=f"{RESULTS_DIR}/summarise_tests_available/summarise_status.csv",
    log:
        f"{LOGS_DIR}/summarise_tests_available.log",
    params:
        src_dir=f"{RESULTS_DIR}/univariate_test",
        module_normality_dir=f"{RESULTS_DIR}/test_normality_modules",
        parametric_method=PARAMETRIC_METHOD,
        nonparametric_method=NONPARAMETRIC_METHOD,
    shell:
        r"""
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/summarise_tests.R \
          --src_dir {params.src_dir:q} \
          --module_normality_dir {params.module_normality_dir:q} \
          --out_dir "$(dirname {output.summary:q})" \
          --parametric_method {params.parametric_method:q} \
          --nonparametric_method {params.nonparametric_method:q} \
          --available_only true \
          > {log:q} 2>&1
        """


rule summarise_tests_members_available:
    output:
        summary=f"{RESULTS_DIR}/summarise_tests_members_available/summary_df.csv",
        status=f"{RESULTS_DIR}/summarise_tests_members_available/summarise_status.csv",
    log:
        f"{LOGS_DIR}/summarise_tests_members_available.log",
    params:
        src_dir=f"{RESULTS_DIR}/univariate_test_members",
        normality_dir=f"{RESULTS_DIR}/test_normality",
        parametric_method=PARAMETRIC_METHOD,
        nonparametric_method=NONPARAMETRIC_METHOD,
    shell:
        r"""
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/summarise_tests.R \
          --src_dir {params.src_dir:q} \
          --module_normality_dir {params.normality_dir:q} \
          --out_dir "$(dirname {output.summary:q})" \
          --parametric_method {params.parametric_method:q} \
          --nonparametric_method {params.nonparametric_method:q} \
          --available_only true \
          > {log:q} 2>&1
        """
