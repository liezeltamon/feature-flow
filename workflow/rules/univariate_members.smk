rule univariate_test_members:
    input:
        feature_path=feature_input,
        sample_metadata_path=metadata_input,
    output:
        result=f"{RESULTS_DIR}/univariate_test_members/{{dataset_id}}/{{method}}/result_matrices.rds",
    log:
        f"{LOGS_DIR}/univariate_test_members/{{dataset_id}}.{{method}}.log",
    params:
        contrast=join_space(config["univariate"]["contrast"]),
        group_key_levels=join_space(
            config["univariate"].get("member_group_order", config["univariate"]["group_order"])
        ),
        random_effects=join_space(config["univariate"]["random_effects"]),
        fixed_effects=join_space(config["univariate"]["fixed_effects"]),
        additional_contrasts_json=json.dumps(
            config["univariate"].get(
                "member_additional_contrasts",
                config["univariate"].get("additional_contrasts", {})
            )
        ),
        equivalence_deltas=join_comma(
            config["univariate"].get("equivalence_deltas", [0.05, 0.1, 0.2])
        ),
        compute_model_assumptions=lambda wildcards: bool_str(
            wildcards.method == config.get("summarise_tests", {}).get("parametric_method", "lmer")
        ),
        random_effects_arg=(
            ""
            if not config["univariate"]["random_effects"]
            else f"--random_effects {join_space(config['univariate']['random_effects'])}"
        ),
        fixed_effects_arg=(
            ""
            if not config["univariate"]["fixed_effects"]
            else f"--fixed_effects {join_space(config['univariate']['fixed_effects'])}"
        ),
    shell:
        r"""
        mkdir -p "$(dirname {output.result:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/univariate_test.R \
          --feature_path {input.feature_path:q} \
          --sample_metadata_path {input.sample_metadata_path:q} \
          --out_dir "$(dirname {output.result:q})" \
          --sample_key {SAMPLE_KEY:q} \
          --contrast {params.contrast} \
          --group_key_levels {params.group_key_levels} \
          {params.random_effects_arg} \
          {params.fixed_effects_arg} \
          --method {wildcards.method:q} \
          --additional_contrasts_json {params.additional_contrasts_json:q} \
          --equivalence_deltas {params.equivalence_deltas:q} \
          --compute_model_assumptions {params.compute_model_assumptions:q} \
          --missing_response_policy drop_feature_sample \
          --min_samples_per_group 3 \
          > {log:q} 2>&1
        """
