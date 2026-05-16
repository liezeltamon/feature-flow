rule plot_features_members:
    input:
        plot_feature_path=plot_feature_input,
        sample_metadata_path=metadata_input,
    output:
        marker=f"{RESULTS_DIR}/plot_features_members/{{dataset_id}}/m1/features.pdf",
    log:
        f"{LOGS_DIR}/plot_features_members/{{dataset_id}}.log",
    params:
        plot_group_order=join_space(
            config.get("plot_features", {}).get("group_order", config["univariate"]["group_order"])
        ),
    shell:
        r"""
        mkdir -p "$(dirname {output.marker:q})/.." "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/plot_features.R \
          --plot_feature_path {input.plot_feature_path:q} \
          --sample_metadata_path {input.sample_metadata_path:q} \
          --out_dir "$(dirname {output.marker:q})/.." \
          --sample_key {SAMPLE_KEY:q} \
          --group_key {GROUP_KEY:q} \
          --module_prefix "" \
          --plot_group_order {params.plot_group_order} \
          > {log:q} 2>&1
        """
