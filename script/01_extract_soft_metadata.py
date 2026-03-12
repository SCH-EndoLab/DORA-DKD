#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOFT_DIR = ROOT / "data" / "raw" / "geo_soft"
OUT_DIR = ROOT / "data" / "metadata" / "soft_samples"
SUMMARY_FILE = OUT_DIR / "sample_sheet_summary.tsv"


def normalize_key(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def parse_soft_file(path: Path) -> list[dict[str, list[str]]]:
    records: list[dict[str, list[str]]] = []
    current: dict[str, list[str]] | None = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = defaultdict(list)
                current["sample_id"].append(line.split(" = ", 1)[1].strip())
                continue

            if current is None or not line.startswith("!Sample_"):
                continue

            payload = line[len("!Sample_") :]
            if " = " not in payload:
                continue

            key, value = payload.split(" = ", 1)
            current[key.strip()].append(value.strip())

    if current is not None:
        records.append(current)

    return records


def fold_characteristics(values: list[str]) -> dict[str, str]:
    folded: dict[str, list[str]] = defaultdict(list)

    for value in values:
        if ":" in value:
            key, item = value.split(":", 1)
            folded[normalize_key(key)].append(item.strip())
        else:
            folded["characteristics_unparsed"].append(value.strip())

    return {key: " | ".join(items) for key, items in folded.items()}


def flatten_record(dataset_id: str, record: dict[str, list[str]]) -> dict[str, str]:
    flat: dict[str, str] = {"dataset_id": dataset_id}

    preferred_keys = [
        "sample_id",
        "title",
        "geo_accession",
        "organism_ch1",
        "source_name_ch1",
        "platform_id",
        "type",
        "relation",
        "supplementary_file",
    ]

    for key in preferred_keys:
        if key in record:
            flat[key] = " | ".join(record[key])

    characteristics = fold_characteristics(record.get("characteristics_ch1", []))
    for key, value in characteristics.items():
        flat[f"char_{key}"] = value

    if "characteristics_ch1" in record:
        flat["characteristics_raw"] = " | ".join(record["characteristics_ch1"])

    return flat


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return

    all_keys: list[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                all_keys.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []

    for soft_path in sorted(SOFT_DIR.glob("*_family.soft.gz")):
        dataset_id = soft_path.name.replace("_family.soft.gz", "")
        records = parse_soft_file(soft_path)
        rows = [flatten_record(dataset_id, record) for record in records]
        out_path = OUT_DIR / f"{dataset_id}_samples.tsv"
        write_tsv(out_path, rows)

        summary_rows.extend(rows)
        print(f"{dataset_id}: parsed {len(rows)} samples -> {out_path}")

    write_tsv(SUMMARY_FILE, summary_rows)
    print(f"summary: wrote {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
