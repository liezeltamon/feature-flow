rule resolved_dataset_config:
    output:
        f"{RESULTS_DIR}/resolved_dataset_config.csv",
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        fieldnames = [
            "dataset_id",
            "min_cluster_size",
            "n_repeats",
            "filter_features_by_subset_unique_values",
            "subset_min_unique_values",
            "matched_rules",
            "has_exact_override",
        ]
        with open(output[0], "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for dataset_id in DATASET_IDS:
                modules = module_config(dataset_id)
                writer.writerow({
                    "dataset_id": dataset_id,
                    "min_cluster_size": modules["min_cluster_size"],
                    "n_repeats": modules["n_repeats"],
                    "filter_features_by_subset_unique_values": modules.get(
                        "filter_features_by_subset_unique_values", False
                    ),
                    "subset_min_unique_values": modules.get("subset_min_unique_values"),
                    "matched_rules": ";".join(matched_dataset_override_rule_names(dataset_id)),
                    "has_exact_override": dataset_id in config.get("dataset_overrides", {}),
                })
