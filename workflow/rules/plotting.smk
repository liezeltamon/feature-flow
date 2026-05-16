rule plot_features:
    input:
        plot_feature_path=plot_feature_input,
        module_feature_path=lambda wildcards: canonical_module_table_path(wildcards.dataset_id),
        eigengenes_path=rules.modules_hc.output.eigengenes,
        sample_metadata_path=plot_metadata_input,
    output:
        marker=f"{RESULTS_DIR}/plot_features/{{dataset_id}}/{SUBSET_NAME}/variance_explained.svg",
    log:
        f"{LOGS_DIR}/plot_features/{{dataset_id}}.{SUBSET_NAME}.log",
    params:
        plot_group_order=join_space(
            config.get("plot_features", {}).get("group_order", config["univariate"]["group_order"])
        ),
        module_prefix=lambda wildcards: f"{SUBSET_NAME}..{wildcards.dataset_id}..m..",
        module_annotations_arg=lambda wildcards: (
            ""
            if canonical_module_annotations_path(wildcards.dataset_id) is None
            else f"--module_annotations_path {canonical_module_annotations_path(wildcards.dataset_id)}"
        ),
    shell:
        r"""
        mkdir -p "$(dirname {output.marker:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/plot_features.R \
          --plot_feature_path {input.plot_feature_path:q} \
          --module_feature_path {input.module_feature_path:q} \
          --eigengenes_path {input.eigengenes_path:q} \
          --sample_metadata_path {input.sample_metadata_path:q} \
          --out_dir "$(dirname {output.marker:q})" \
          --sample_key {SAMPLE_KEY:q} \
          --group_key {GROUP_KEY:q} \
          --module_prefix {params.module_prefix:q} \
          {params.module_annotations_arg} \
          --plot_group_order {params.plot_group_order} \
          > {log:q} 2>&1
        """
