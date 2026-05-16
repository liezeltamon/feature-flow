rule test_normality_modules:
    input:
        feature_path=lambda wildcards: canonical_module_table_path(wildcards.dataset_id),
    output:
        results=f"{RESULTS_DIR}/test_normality_modules/{{dataset_id}}/{SUBSET_NAME}/normality_test_results.csv",
        summary=f"{RESULTS_DIR}/test_normality_modules/{{dataset_id}}/{SUBSET_NAME}/normality_summary.csv",
    log:
        f"{LOGS_DIR}/test_normality_modules/{{dataset_id}}.{SUBSET_NAME}.log",
    params:
        alpha=config["normality"]["alpha"],
        correct=config["normality"]["correct"],
        test=config["normality"]["test"],
    shell:
        r"""
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/test_normality_wrapper.R \
          --feature_path {input.feature_path:q} \
          --sample_key {SAMPLE_KEY:q} \
          --out_dir "$(dirname {output.summary:q})" \
          --alpha {params.alpha} \
          --correct {params.correct:q} \
          --test {params.test:q} \
          > {log:q} 2>&1
        """
