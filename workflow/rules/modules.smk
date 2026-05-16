if MODULES_HC_SLURM_PARTITION is None:
    rule modules_hc:
        input:
            feature_path=rules.preprocess_impute.output.processed,
            preprocess_summary=rules.preprocess_impute.output.summary,
            metadata_path=metadata_input,
        output:
            module_table=f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/bulk_x_features.csv",
            eigengenes=f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/module_eigengenes_hc.rds",
            module_members=directory(f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/module_members"),
        threads: RULE_THREADS.get("modules_hc", 1)
        resources:
            mem_mb=rule_resource("modules_hc", "mem_mb", 16000),
        log:
            f"{LOGS_DIR}/modules_hc/{{dataset_id}}.{SUBSET_NAME}.log",
        params:
            min_cluster_size=lambda wildcards: module_param(wildcards.dataset_id, "min_cluster_size"),
            n_repeats=lambda wildcards: module_param(wildcards.dataset_id, "n_repeats"),
            filter_features_by_subset_unique_values=lambda wildcards: bool_str(
                optional_module_param(wildcards.dataset_id, "filter_features_by_subset_unique_values") or False
            ),
            subset_key_arg=lambda wildcards: subset_key_arg(wildcards.dataset_id),
            subset_values_arg=lambda wildcards: subset_values_arg(wildcards.dataset_id),
            subset_min_unique_values_arg=lambda wildcards: subset_min_unique_values_arg(wildcards.dataset_id),
            resolved_module_config_message=lambda wildcards: resolved_module_config_message(wildcards.dataset_id),
        shell:
            r"""
            mkdir -p "$(dirname {output.module_table:q})" "$(dirname {log:q})"
            {{
              echo {params.resolved_module_config_message:q}
              {RSCRIPT_BIN:q} code/feature_module_reduction/modules_hc.R \
                --src_id {wildcards.dataset_id:q} \
                --feature_path {input.feature_path:q} \
                --preprocess_summary_path {input.preprocess_summary:q} \
                --metadata_path {input.metadata_path:q} \
                --out_dir {MODULES_ROOT:q} \
                --n_cores {threads} \
                --sample_key {SAMPLE_KEY:q} \
                {params.subset_key_arg} \
                {params.subset_values_arg} \
                --subset_name {SUBSET_NAME:q} \
                --min_cluster_size {params.min_cluster_size} \
                --n_repeats {params.n_repeats} \
                --filter_features_by_subset_unique_values {params.filter_features_by_subset_unique_values:q} \
                {params.subset_min_unique_values_arg}
            }} > {log:q} 2>&1
            """
else:
    rule modules_hc:
        input:
            feature_path=rules.preprocess_impute.output.processed,
            preprocess_summary=rules.preprocess_impute.output.summary,
            metadata_path=metadata_input,
        output:
            module_table=f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/bulk_x_features.csv",
            eigengenes=f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/module_eigengenes_hc.rds",
            module_members=directory(f"{RESULTS_DIR}/modules_hc/{{dataset_id}}/{SUBSET_NAME}/module_members"),
        threads: RULE_THREADS.get("modules_hc", 1)
        resources:
            mem_mb=rule_resource("modules_hc", "mem_mb", 16000),
            slurm_partition=MODULES_HC_SLURM_PARTITION,
        log:
            f"{LOGS_DIR}/modules_hc/{{dataset_id}}.{SUBSET_NAME}.log",
        params:
            min_cluster_size=lambda wildcards: module_param(wildcards.dataset_id, "min_cluster_size"),
            n_repeats=lambda wildcards: module_param(wildcards.dataset_id, "n_repeats"),
            filter_features_by_subset_unique_values=lambda wildcards: bool_str(
                optional_module_param(wildcards.dataset_id, "filter_features_by_subset_unique_values") or False
            ),
            subset_key_arg=lambda wildcards: subset_key_arg(wildcards.dataset_id),
            subset_values_arg=lambda wildcards: subset_values_arg(wildcards.dataset_id),
            subset_min_unique_values_arg=lambda wildcards: subset_min_unique_values_arg(wildcards.dataset_id),
            resolved_module_config_message=lambda wildcards: resolved_module_config_message(wildcards.dataset_id),
        shell:
            r"""
            mkdir -p "$(dirname {output.module_table:q})" "$(dirname {log:q})"
            {{
              echo {params.resolved_module_config_message:q}
              {RSCRIPT_BIN:q} code/feature_module_reduction/modules_hc.R \
                --src_id {wildcards.dataset_id:q} \
                --feature_path {input.feature_path:q} \
                --preprocess_summary_path {input.preprocess_summary:q} \
                --metadata_path {input.metadata_path:q} \
                --out_dir {MODULES_ROOT:q} \
                --n_cores {threads} \
                --sample_key {SAMPLE_KEY:q} \
                {params.subset_key_arg} \
                {params.subset_values_arg} \
                --subset_name {SUBSET_NAME:q} \
                --min_cluster_size {params.min_cluster_size} \
                --n_repeats {params.n_repeats} \
                --filter_features_by_subset_unique_values {params.filter_features_by_subset_unique_values:q} \
                {params.subset_min_unique_values_arg}
            }} > {log:q} 2>&1
            """


rule meaning_module_genes:
    input:
        module_feature_path=rules.modules_hc.output.module_table,
        modules_dir=rules.modules_hc.output.module_members,
    output:
        annotated_table=f"{RESULTS_DIR}/meaning_module_genes/{{dataset_id}}/{SUBSET_NAME}/bulk_x_features.csv",
        annotations=f"{RESULTS_DIR}/meaning_module_genes/{{dataset_id}}/{SUBSET_NAME}/module_annotations.csv",
    log:
        f"{LOGS_DIR}/meaning_module_genes/{{dataset_id}}.{SUBSET_NAME}.log",
    params:
        backend=config.get("module_meaning", {}).get("backend", "gprofiler"),
        gmt_path=config.get("module_meaning", {}).get("gmt_path", ""),
        min_hits=config.get("module_meaning", {}).get("min_hits", 2),
        max_frac=config.get("module_meaning", {}).get("max_frac", 0.5),
        organism=config["module_meaning"]["organism"],
        top_n=config["module_meaning"]["top_n"],
        use_only_important_loading=bool_str(config["module_meaning"]["use_only_important_loading"]),
        feature_id_separator=MODULE_MEMBER_FEATURE_ID_SEPARATOR,
        duplicate_enrichment_id_policy=MODULE_MEMBER_DUPLICATE_POLICY,
    shell:
        r"""
        mkdir -p "$(dirname {output.annotated_table:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/meaning_module_genes.R \
          --module_feature_path {input.module_feature_path:q} \
          --modules_dir {input.modules_dir:q} \
          --out_dir "$(dirname {output.annotated_table:q})" \
          --backend {params.backend:q} \
          --gmt_path {params.gmt_path:q} \
          --min_hits {params.min_hits} \
          --max_frac {params.max_frac} \
          --organism {params.organism:q} \
          --top_n {params.top_n} \
          --use_only_important_loading {params.use_only_important_loading:q} \
          --feature_id_separator {params.feature_id_separator:q} \
          --duplicate_enrichment_id_policy {params.duplicate_enrichment_id_policy:q} \
          > {log:q} 2>&1
        """


rule meaning_module_llm:
    input:
        module_feature_path=rules.modules_hc.output.module_table,
        modules_dir=rules.modules_hc.output.module_members,
    output:
        annotated_table=f"{RESULTS_DIR}/meaning_module_llm/{{dataset_id}}/{SUBSET_NAME}/bulk_x_features.csv",
        annotations=f"{RESULTS_DIR}/meaning_module_llm/{{dataset_id}}/{SUBSET_NAME}/module_annotations.csv",
    log:
        f"{LOGS_DIR}/meaning_module_llm/{{dataset_id}}.{SUBSET_NAME}.log",
    params:
        top_n=config.get("module_meaning_llm", {}).get("top_n", 5),
        model=config.get("module_meaning_llm", {}).get("model", "gpt-5-nano"),
        use_only_important_loading_arg=(
            "--use_only_important_loading"
            if config.get("module_meaning_llm", {}).get("use_only_important_loading", False)
            else ""
        ),
        feature_id_separator=MODULE_MEMBER_FEATURE_ID_SEPARATOR,
        duplicate_enrichment_id_policy=MODULE_MEMBER_DUPLICATE_POLICY,
    shell:
        r"""
        mkdir -p "$(dirname {output.annotated_table:q})" "$(dirname {log:q})"
        {PYTHON_BIN:q} code/feature_module_reduction/meaning_module_llm.py \
          --module_feature_path {input.module_feature_path:q} \
          --modules_dir {input.modules_dir:q} \
          --out_dir "$(dirname {output.annotated_table:q})" \
          --top_n {params.top_n} \
          --model {params.model:q} \
          --feature_id_separator {params.feature_id_separator:q} \
          --duplicate_enrichment_id_policy {params.duplicate_enrichment_id_policy:q} \
          {params.use_only_important_loading_arg} \
          > {log:q} 2>&1
        """


rule module_set_enrichment:
    input:
        modules_dir=rules.modules_hc.output.module_members,
        gmt_path=lambda wildcards: MODULE_SET_ENRICHMENT["gmt_path"],
    output:
        enrichment=f"{RESULTS_DIR}/module_set_enrichment/{{dataset_id}}/{SUBSET_NAME}/mhg_enrichment.csv",
    log:
        f"{LOGS_DIR}/module_set_enrichment/{{dataset_id}}.{SUBSET_NAME}.log",
    params:
        feature_id_format=MODULE_SET_ENRICHMENT.get("feature_id_format", "module_member_feature_ids"),
        min_hits=MODULE_SET_ENRICHMENT.get("min_hits", 2),
        max_frac=MODULE_SET_ENRICHMENT.get("max_frac", 0.5),
        use_only_important_loading=bool_str(MODULE_SET_ENRICHMENT.get("use_only_important_loading", True)),
        feature_id_separator=MODULE_MEMBER_FEATURE_ID_SEPARATOR,
        duplicate_enrichment_id_policy=MODULE_MEMBER_DUPLICATE_POLICY,
    shell:
        r"""
        mkdir -p "$(dirname {output.enrichment:q})" "$(dirname {log:q})"
        {RSCRIPT_BIN:q} code/feature_module_reduction/module_set_enrichment.R \
          --modules_dir {input.modules_dir:q} \
          --gmt_path {input.gmt_path:q} \
          --out_dir "$(dirname {output.enrichment:q})" \
          --feature_id_format {params.feature_id_format:q} \
          --feature_id_separator {params.feature_id_separator:q} \
          --duplicate_enrichment_id_policy {params.duplicate_enrichment_id_policy:q} \
          --min_hits {params.min_hits} \
          --max_frac {params.max_frac} \
          --use_only_important_loading {params.use_only_important_loading:q} \
          > {log:q} 2>&1
        """
