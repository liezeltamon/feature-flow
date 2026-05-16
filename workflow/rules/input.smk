if INPUT_MODE == "combine":
    rule combine_inputs:
        params:
            input_dir=config["input_root"],
            combine_mode=COMBINE_MODE,
            feature_filename=FEATURE_FILENAME,
            metadata_filename=METADATA_FILENAME,
        output:
            feature_path=COMBINE_FEATURE_PATH,
            metadata_path=COMBINE_METADATA_PATH,
            manifest=f"{COMBINE_ROOT}/input_manifest.csv",
        log:
            f"{LOGS_DIR}/combine_inputs/{COMBINE_ID}.log",
        shell:
            r"""
            mkdir -p "$(dirname {output.feature_path:q})" "$(dirname {log:q})"
            {RSCRIPT_BIN:q} code/feature_module_reduction/combine_tables.R \
              --mode {params.combine_mode:q} \
              --input_dir {params.input_dir:q} \
              --feature_filename {params.feature_filename:q} \
              --metadata_filename {params.metadata_filename:q} \
              --sample_key {SAMPLE_KEY:q} \
              --out_dir {COMBINE_ROOT:q} \
              > {log:q} 2>&1
            """
