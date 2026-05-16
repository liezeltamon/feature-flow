shell.executable("/bin/bash")

import os
import json
import csv
import re


def join_space(values):
    if values is None:
        return ""
    if isinstance(values, str):
        values = [values]
    values = [str(v) for v in values if v is not None and str(v) != ""]
    return " ".join(values)


def join_comma(values):
    if values is None:
        return ""
    if isinstance(values, str):
        return values
    return ",".join(str(v) for v in values)


def bool_str(value):
    return "true" if value else "false"


def deep_merge(base, override):
    result = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


RESULTS_DIR = config["results_dir"]
LOGS_DIR = config["logs_dir"]
RSCRIPT_BIN = config["rscript_bin"]
PYTHON_BIN = config.get("python_bin", "python3")
SAMPLE_KEY = config["sample_key"]
GROUP_KEY = config["group_key"]
INDIVIDUAL_KEY = config["individual_key"]
SUBSET_NAME = config["modules"]["subset_name"]
RULE_THREADS = config.get("threads", {})
RULE_RESOURCES = config.get("resources", {})
MODULE_MEANING_GENES_ENABLED = config["module_meaning"]["enabled"]
MODULE_MEANING_LLM_ENABLED = config.get("module_meaning_llm", {}).get("enabled", False)
MODULE_MEMBER_FEATURE_IDS = config.get("module_member_feature_ids", {})
MODULE_MEMBER_FEATURE_ID_SEPARATOR = MODULE_MEMBER_FEATURE_IDS.get("id_separator", "__")
MODULE_MEMBER_DUPLICATE_POLICY = MODULE_MEMBER_FEATURE_IDS.get(
    "duplicate_enrichment_id_policy",
    "max_abs_loading",
)
MODULE_SET_ENRICHMENT = config.get("module_set_enrichment", {})
MODULE_SET_ENRICHMENT_ENABLED = MODULE_SET_ENRICHMENT.get("enabled", False)
INPUT_MODE = config.get("input_mode", "single")
COMBINE_MODE = config.get("combine_mode", "row_bind")
FEATURE_FILENAME = config.get("feature_filename", "feature_table.csv")
METADATA_FILENAME = config.get("metadata_filename", "sample_metadata.csv")
PLOT_INPUT_ROOT = config.get("plot_input_root", None)
PLOT_FEATURE_FILENAME = config.get("plot_feature_filename", FEATURE_FILENAME)
PLOT_METADATA_FILENAME = config.get("plot_metadata_filename", METADATA_FILENAME)
COMBINE_ID = config.get("combine", {}).get("combine_id", "combined")
GENE_DATASET_MARKERS = ("var_genes",)

if INPUT_MODE == "combine" and COMBINE_MODE not in ("row_bind", "column_bind"):
    raise ValueError(f"Unsupported combine_mode: {COMBINE_MODE}")


def rule_resource(rule_name, resource_name, default):
    return RULE_RESOURCES.get(rule_name, {}).get(resource_name, default)


def optional_rule_resource(rule_name, resource_name):
    value = RULE_RESOURCES.get(rule_name, {}).get(resource_name)
    if value in (None, ""):
        return None
    return value


def discover_dataset_ids():
    if INPUT_MODE == "single":
        return [config["src_id"]]

    if INPUT_MODE == "combine":
        return [COMBINE_ID]

    if INPUT_MODE == "batch":
        dataset_ids = []
        for name in sorted(os.listdir(config["input_root"])):
            path = os.path.join(config["input_root"], name)
            if (
                os.path.isdir(path)
                and os.path.exists(os.path.join(path, FEATURE_FILENAME))
                and os.path.exists(os.path.join(path, METADATA_FILENAME))
            ):
                dataset_ids.append(name)
        if len(dataset_ids) == 0:
            raise ValueError("No valid dataset subdirectories found for batch mode")
        return dataset_ids

    raise ValueError(f"Unsupported input_mode: {INPUT_MODE}")


DATASET_IDS = discover_dataset_ids()


def matching_dataset_override_rules(dataset_id):
    matched = []
    for rule in config.get("dataset_override_rules", []):
        pattern = rule.get("match")
        if pattern is not None and re.search(pattern, dataset_id):
            matched.append(rule)
    return matched


def resolve_dataset_config(dataset_id):
    resolved = dict(config)
    for rule in matching_dataset_override_rules(dataset_id):
        override = {k: v for k, v in rule.items() if k not in ("name", "match")}
        resolved = deep_merge(resolved, override)
    exact_override = config.get("dataset_overrides", {}).get(dataset_id, {})
    if exact_override:
        resolved = deep_merge(resolved, exact_override)
    return resolved


def matched_dataset_override_rule_names(dataset_id):
    names = []
    for i, rule in enumerate(matching_dataset_override_rules(dataset_id), start=1):
        names.append(str(rule.get("name", f"rule_{i}")))
    return names


def module_config(dataset_id):
    return resolve_dataset_config(dataset_id)["modules"]


def module_param(dataset_id, key):
    return module_config(dataset_id)[key]


def optional_module_param(dataset_id, key):
    return module_config(dataset_id).get(key)


def subset_key_arg(dataset_id):
    subset_key = optional_module_param(dataset_id, "subset_key")
    if subset_key is None:
        return ""
    return f"--subset_key {subset_key}"


def subset_values_arg(dataset_id):
    subset_values = optional_module_param(dataset_id, "subset_values")
    if not subset_values:
        return ""
    return f"--subset_values {join_comma(subset_values)}"


def subset_min_unique_values_arg(dataset_id):
    subset_min_unique_values = optional_module_param(dataset_id, "subset_min_unique_values")
    if subset_min_unique_values is None:
        return ""
    return f"--subset_min_unique_values {subset_min_unique_values}"


def resolved_module_config_message(dataset_id):
    modules = module_config(dataset_id)
    matched_rules = ";".join(matched_dataset_override_rule_names(dataset_id))
    return (
        f"Resolved module config for {dataset_id}: "
        f"min_cluster_size={modules['min_cluster_size']}; "
        f"n_repeats={modules['n_repeats']}; "
        "filter_features_by_subset_unique_values="
        f"{modules.get('filter_features_by_subset_unique_values', False)}; "
        f"subset_min_unique_values={modules.get('subset_min_unique_values')}; "
        f"matched_rules={matched_rules}; "
        f"has_exact_override={dataset_id in config.get('dataset_overrides', {})}"
    )


def is_gene_like_dataset(dataset_id):
    return any(marker in dataset_id for marker in GENE_DATASET_MARKERS)


GENE_LIKE_DATASET_IDS = [dataset_id for dataset_id in DATASET_IDS if is_gene_like_dataset(dataset_id)]

COMBINE_ROOT = f"{RESULTS_DIR}/combine_inputs/{COMBINE_ID}"
COMBINE_FEATURE_PATH = f"{COMBINE_ROOT}/feature_table.csv"
COMBINE_METADATA_PATH = f"{COMBINE_ROOT}/sample_metadata.csv"


def canonical_module_table_path(dataset_id):
    if MODULE_MEANING_GENES_ENABLED and is_gene_like_dataset(dataset_id):
        return f"{RESULTS_DIR}/meaning_module_genes/{dataset_id}/{SUBSET_NAME}/bulk_x_features.csv"
    if MODULE_MEANING_LLM_ENABLED:
        return f"{RESULTS_DIR}/meaning_module_llm/{dataset_id}/{SUBSET_NAME}/bulk_x_features.csv"
    return f"{RESULTS_DIR}/modules_hc/{dataset_id}/{SUBSET_NAME}/bulk_x_features.csv"


def canonical_module_annotations_path(dataset_id):
    if MODULE_MEANING_GENES_ENABLED and is_gene_like_dataset(dataset_id):
        return f"{RESULTS_DIR}/meaning_module_genes/{dataset_id}/{SUBSET_NAME}/module_annotations.csv"
    if MODULE_MEANING_LLM_ENABLED:
        return f"{RESULTS_DIR}/meaning_module_llm/{dataset_id}/{SUBSET_NAME}/module_annotations.csv"
    return None


def feature_input_path(dataset_id):
    if INPUT_MODE == "single":
        return config["feature_path"]
    if INPUT_MODE == "combine":
        return COMBINE_FEATURE_PATH
    return os.path.join(config["input_root"], dataset_id, FEATURE_FILENAME)


def metadata_input_path(dataset_id):
    if INPUT_MODE == "single":
        return config["sample_metadata_path"]
    if INPUT_MODE == "combine":
        return COMBINE_METADATA_PATH
    return os.path.join(config["input_root"], dataset_id, METADATA_FILENAME)


def feature_input(wildcards):
    return feature_input_path(wildcards.dataset_id)


def metadata_input(wildcards):
    return metadata_input_path(wildcards.dataset_id)


def plot_feature_input_path(dataset_id):
    if PLOT_INPUT_ROOT is not None:
        return os.path.join(PLOT_INPUT_ROOT, dataset_id, PLOT_FEATURE_FILENAME)
    return feature_input_path(dataset_id)


def plot_metadata_input_path(dataset_id):
    if PLOT_INPUT_ROOT is not None:
        return os.path.join(PLOT_INPUT_ROOT, dataset_id, PLOT_METADATA_FILENAME)
    return metadata_input_path(dataset_id)


def plot_feature_input(wildcards):
    return plot_feature_input_path(wildcards.dataset_id)


def plot_metadata_input(wildcards):
    return plot_metadata_input_path(wildcards.dataset_id)


MODULES_ROOT = f"{RESULTS_DIR}/modules_hc"
EVALUATE_PREPROCESSING_ROOT = f"{RESULTS_DIR}/evaluate_preprocessing/{config['preprocess']['error_metric']}"
MODULES_HC_SLURM_PARTITION = optional_rule_resource("modules_hc", "slurm_partition")

if INPUT_MODE == "batch":
    EVALUATE_FEATURE_TABLE_DIR = config["input_root"]
    EVALUATE_FEATURE_FILENAME = FEATURE_FILENAME
elif INPUT_MODE == "single":
    EVALUATE_FEATURE_TABLE_DIR = os.path.dirname(config["feature_path"])
    EVALUATE_FEATURE_FILENAME = os.path.basename(config["feature_path"])
elif INPUT_MODE == "combine":
    EVALUATE_FEATURE_TABLE_DIR = os.path.dirname(COMBINE_FEATURE_PATH)
    EVALUATE_FEATURE_FILENAME = os.path.basename(COMBINE_FEATURE_PATH)
else:
    raise ValueError(f"Unsupported input_mode: {INPUT_MODE}")
