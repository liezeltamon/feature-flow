#!/usr/bin/env python3

# %% Setup

import argparse
import json
import re
from pathlib import Path

import pandas as pd
from openai import OpenAI


# %% Parameters

parser = argparse.ArgumentParser(description="Add LLM-based meaning labels to module columns")
parser.add_argument("--module_feature_path", required=True)
parser.add_argument("--modules_dir", required=True)
parser.add_argument("--out_dir", required=True)
parser.add_argument("--top_n", type=int, default=5)
parser.add_argument("--model", default="gpt-5-nano")
parser.add_argument("--use_only_important_loading", action="store_true")
parser.add_argument("--feature_id_separator", default="__")
parser.add_argument("--duplicate_enrichment_id_policy", default="max_abs_loading")
args = parser.parse_args()

if not args.feature_id_separator:
    raise ValueError("feature_id_separator must be a non-empty string")

if args.duplicate_enrichment_id_policy != "max_abs_loading":
    raise ValueError(
        "duplicate_enrichment_id_policy must be max_abs_loading; got: "
        f"{args.duplicate_enrichment_id_policy}"
    )


# %% Helpers

def sanitize_label(x):
    x = str(x).replace(" ", "_")
    x = re.sub(r"[^A-Za-z0-9_.-]", "_", x)
    x = re.sub(r"_+", "_", x)
    return x.strip("_")


def sanitize_id_part(x):
    x = str(x).replace(" ", "_")
    x = re.sub(r"[^A-Za-z0-9_.-]", "_", x)
    return x.strip("_")


def extract_feature_id(feature, feature_id_separator):
    return str(feature).split(feature_id_separator, maxsplit=1)[0]


def parse_module_id(module_id):
    parts = str(module_id).split("..")
    if len(parts) < 4 or parts[-2] != "m" or not parts[-1]:
        raise ValueError(
            "Expected module_id format <subset>..<dataset_id>..m..<number>; "
            f"got: {module_id}"
        )

    dataset_id = "..".join(parts[1:-2])
    if not dataset_id:
        raise ValueError(f"Could not parse dataset_id from module_id: {module_id}")

    return dataset_id, f"m..{parts[-1]}"


def build_renamed_module(module_id, meaning_label):
    dataset_id, module_short_id = parse_module_id(module_id)
    return "__".join([
        sanitize_label(meaning_label),
        sanitize_id_part(dataset_id),
        sanitize_id_part(module_short_id),
    ])


def collect_fallback_feature_ids(df, feature_id_separator):
    fallback_df = df.copy()
    fallback_df["feature_id"] = fallback_df["feature"].map(
        lambda x: extract_feature_id(x, feature_id_separator)
    )

    if fallback_df["feature_id"].isna().any() or (fallback_df["feature_id"] == "").any():
        raise ValueError("Blank extracted feature IDs found while building fallback label")

    fallback_df["abs_loading"] = fallback_df["loading"].abs()
    fallback_df = (
        fallback_df.sort_values("abs_loading", ascending=False)
        .drop_duplicates("feature_id", keep="first")
        .sort_values("loading", ascending=False)
    )
    return fallback_df["feature_id"].astype(str).tolist()


def build_prompt(module_id, top_df):
    evidence_lines = [
        f"- feature={row.feature} | loading={row.loading:.4f}"
        for row in top_df.itertuples(index=False)
    ]

    return "\n".join([
        "You are naming a biological feature module.",
        "",
        "Context:",
        "- These features are module members from a refactored workflow based on WGCNA.",
        "- The loading values are PC/module loadings.",
        "- Higher absolute loading means a feature contributes more strongly to the module.",
        "- Use the loading values as weights when deciding on the name.",
        "- Balance specificity with coverage of the top-loading module members. The label should be specific enough to attach biological meaning to this module, while still reflecting the shared signal across the top-loading members as much as possible.",
        "- The name should summarize the shared biological theme of the top-loading module members, with higher-loading members carrying more weight.",
        "- Do not simply list feature names.",
        "- Do not let a low-loading outlier dominate the name.",
        "- Use only the evidence below.",
        "- If the top-loading members do not support a clear shared theme, say that clearly.",
        "",
        "Return valid JSON with keys:",
        "short_label, summary, confidence, ambiguity_note",
        "",
        "Rules for the response:",
        "- short_label: at most 6 words",
        "- summary: one sentence",
        "- confidence: integer from 1 to 5",
        "- ambiguity_note: short phrase, or 'none'",
        "",
        f"Module ID: {module_id}",
        "Evidence:",
        *evidence_lines,
    ])


def call_llm(client, prompt, model):
    response = client.responses.create(
        model=model,
        input=prompt,
        text={"format": {"type": "json_object"}},
    )

    if not response.output_text:
        raise ValueError("No text returned from model")

    return json.loads(response.output_text)


# %% Load inputs

client = OpenAI()

module_feature_df = pd.read_csv(args.module_feature_path)
modules_dir = Path(args.modules_dir)
out_dir = Path(args.out_dir)
modules_out_dir = out_dir / "modules"

out_dir.mkdir(parents=True, exist_ok=True)
modules_out_dir.mkdir(parents=True, exist_ok=True)

module_paths = list(sorted(modules_dir.glob("*.csv")))
module_paths.extend(sorted(modules_dir.glob("**/members.csv")))
module_paths = list(dict.fromkeys(module_paths))

if not module_paths:
    raise ValueError(f"No module CSV files found in {modules_dir}")


# %% Name modules

annotation_rows = []

for module_path in module_paths:
    module_id = module_path.parent.name if module_path.name == "members.csv" else module_path.stem
    print(f"Naming module: {module_id}")

    if module_id not in module_feature_df.columns:
        raise ValueError(f"{module_id} not found in module feature table")

    df = pd.read_csv(module_path).sort_values("loading", ascending=False)

    if "feature" not in df.columns or "loading" not in df.columns:
        raise ValueError(f"{module_path} must contain feature and loading columns")

    if args.use_only_important_loading:
        if "important_loading" not in df.columns:
            raise ValueError(f"{module_path} missing important_loading column")
        keep = df["important_loading"]
        if keep.dtype != bool:
            keep = keep.astype(str).str.lower().isin(["true", "t", "1"])
        df = df.loc[keep].copy()

    if df.empty:
        raise ValueError(f"No module members left after filtering for {module_id}")

    top_df = df.head(args.top_n).copy()

    module_out_dir = modules_out_dir / module_id
    module_out_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(module_out_dir / "ranked_features.csv", index=False)
    top_df.to_csv(module_out_dir / "top_features.csv", index=False)

    # LLM receives all genes (full df, already filtered by use_only_important_loading)
    prompt = build_prompt(module_id, df)
    (module_out_dir / "prompt.txt").write_text(prompt)

    # top_n only used for the fallback label; sorted alphabetically for consistency
    fallback_feature_ids = collect_fallback_feature_ids(
        df,
        args.feature_id_separator,
    )
    fallback_feature_ids = sorted(fallback_feature_ids[: args.top_n])
    fallback_label = sanitize_label("_".join(fallback_feature_ids))
    short_label = fallback_label
    summary = "LLM naming failed; using fallback from top features."
    confidence = None
    ambiguity_note = "api_error"

    try:
        result = call_llm(client, prompt, args.model)
        (module_out_dir / "response.json").write_text(json.dumps(result, indent=2))

        short_label = result.get("short_label") or short_label
        summary = result.get("summary", summary)
        confidence = result.get("confidence", confidence)
        ambiguity_note = result.get("ambiguity_note", "none")
    except Exception as e:
        (module_out_dir / "response.json").write_text(
            json.dumps({"error": str(e)}, indent=2)
        )

    meaning_label = sanitize_label(short_label) or fallback_label
    renamed_module = build_renamed_module(module_id, meaning_label)

    (module_out_dir / "meaning.txt").write_text(f"{meaning_label}\n")

    annotation_rows.append({
        "module_id": module_id,
        "module_csv_path": str(module_path),
        "top_features": "__".join(top_df["feature"].astype(str).tolist()),
        "top_loadings": "__".join(top_df["loading"].map(lambda x: f"{x:.4f}").tolist()),
        "meaning_label": meaning_label,
        "llm_short_label": short_label,
        "llm_summary": summary,
        "llm_confidence": confidence,
        "llm_ambiguity_note": ambiguity_note,
        "model": args.model,
        "renamed_module": renamed_module,
    })


# %% Write outputs

annotation_df = pd.DataFrame(annotation_rows)

duplicated_renamed_modules = annotation_df.loc[
    annotation_df["renamed_module"].duplicated(),
    "renamed_module",
].drop_duplicates()

if not duplicated_renamed_modules.empty:
    examples = ", ".join(duplicated_renamed_modules.head(5).tolist())
    raise ValueError(f"Duplicated renamed module names found. Examples: {examples}")

annotation_df.to_csv(out_dir / "module_annotations.csv", index=False)

rename_lookup = dict(zip(annotation_df["module_id"], annotation_df["renamed_module"]))
renamed_df = module_feature_df.rename(columns=rename_lookup)
renamed_df.to_csv(out_dir / "bulk_x_features.csv", index=False)
