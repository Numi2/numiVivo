#!/usr/bin/env python3
"""Build the frozen, candidate-blind phase-2 target catalog from three real FASTAs.

Normalization, accession parsing, exact full-sequence deduplication, representative
ordering, and output layout follow the phase-1 catalog builder. This program never
enumerates or reports peptide identities. Comet generates decoys at search time.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO


WORKSPACE = Path(__file__).resolve().parent
LOCK_PATH = WORKSPACE / "SOURCE_LOCK.json"
EXPECTED_LOCK_SHA256 = "b16ddd17d3c819fe7042575d7e990a663712b8f73f90faf410882679129790c2"
SOURCE_ORDER = ("CON", "UPR", "TRC")
SOURCE_RANK = {name: rank for rank, name in enumerate(SOURCE_ORDER, start=1)}
REFERENCE_BUILDER_SHA256 = "303f1e6cebb59a5426f1aaaa9bb3176acfbabb3d7b91b7110c2fa8a1448757d0"
COMPONENT_FILENAMES = {
    "CON": "component-CON-contaminants.fasta",
    "UPR": "component-UPR-uniprot-canonical.fasta",
    "TRC": "component-TRC-transcode-phase2-primary-ge8aa.fasta",
}
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
ACCESSION_SAFE = re.compile(r"^[^\s>]+$")


@dataclass
class Record:
    source_class: str
    source_accession: str
    original_header: str
    sequence: str
    selection_reason: str
    sequence_sha256: str = field(init=False)

    def __post_init__(self) -> None:
        self.sequence_sha256 = hashlib.sha256(self.sequence.encode("ascii")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: Optional[str] = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                if not header:
                    raise ValueError(f"empty FASTA header in {path}:{line_number}")
                chunks = []
            else:
                if header is None:
                    if line.strip():
                        raise ValueError(f"sequence before first header in {path}:{line_number}")
                    continue
                chunks.append("".join(line.split()))
    if header is not None:
        yield header, "".join(chunks).upper()


def source_accession(header: str) -> str:
    token = header.split(None, 1)[0]
    parts = token.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        accession = parts[1]
    else:
        accession = parts[0]
    if not accession or not ACCESSION_SAFE.fullmatch(accession):
        raise ValueError(f"unsafe or empty source accession derived from header: {header!r}")
    return accession


def validate_record(record: Record, origin: Path) -> None:
    if not record.sequence:
        raise ValueError(f"empty sequence for {record.source_accession} in {origin}")
    bad = set(record.sequence) - ALLOWED_AA
    if bad:
        raise ValueError(
            f"invalid amino-acid alphabet for {record.source_accession} in {origin}: {sorted(bad)}"
        )
    if "\t" in record.original_header or "\n" in record.original_header:
        raise ValueError(f"TSV-unsafe original header in {origin}: {record.original_header!r}")


def read_component(
    path: Path, source_class: str, selection_reason: str, minimum_length: int
) -> tuple[list[Record], dict[str, int]]:
    records: list[Record] = []
    raw_count = 0
    excluded_short = 0
    seen_accessions: dict[str, str] = {}
    for header, sequence in iter_fasta(path):
        raw_count += 1
        accession = source_accession(header)
        if accession in seen_accessions:
            raise ValueError(f"duplicate accession {accession!r} in {path}")
        seen_accessions[accession] = sequence
        if len(sequence) < minimum_length:
            excluded_short += 1
            continue
        record = Record(source_class, accession, header, sequence, selection_reason)
        validate_record(record, path)
        records.append(record)
    return records, {"raw_records": raw_count, "excluded_below_min_length": excluded_short}


def raw_substring_bound(sequence_length: int, minimum: int = 8, maximum: int = 15) -> int:
    return sum(max(0, sequence_length - peptide_length + 1) for peptide_length in range(minimum, maximum + 1))


def source_summary(records: list[Record]) -> dict[str, object]:
    residues = Counter()
    for record in records:
        residues.update(record.sequence)
    return {
        "selected_records": len(records),
        "selected_residues": sum(len(record.sequence) for record in records),
        "sequence_alphabet": "".join(sorted(residues)),
        "nonstandard_residue_counts": {
            residue: residues[residue] for residue in sorted(set(residues) - STANDARD_AA)
        },
        "raw_8_15_substring_occurrence_upper_bound": sum(
            raw_substring_bound(len(record.sequence)) for record in records
        ),
    }


def write_fasta(path: Path, records: Iterable[tuple[str, str, str]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for accession, description, sequence in records:
            handle.write(f">{accession}")
            if description:
                handle.write(f" {description}")
            handle.write("\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_lock() -> tuple[dict[str, object], str]:
    lock_sha = sha256_file(LOCK_PATH)
    if lock_sha != EXPECTED_LOCK_SHA256:
        raise ValueError(f"source lock SHA-256 changed: {lock_sha} != {EXPECTED_LOCK_SHA256}")
    lock = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
    if lock["schema"] != "numivivo.melanoma.phase2.candidate_blind_catalog_source_lock.v1":
        raise ValueError("unexpected source-lock schema")
    if tuple(lock["source_order"]) != SOURCE_ORDER or set(lock["sources"]) != set(SOURCE_ORDER):
        raise ValueError("source lock must contain exactly CON, UPR, TRC in frozen order")
    if lock["excluded_sources"] != ["NUO", "SMO"]:
        raise ValueError("source exclusions differ from frozen phase-2 rule")
    policy = lock["search_contract"]
    if policy["catalog_contains_targets_only"] is not True or policy["decoy_search"] != 1:
        raise ValueError("target-only / Comet generated-decoy rule changed")
    if policy["decoy_prefix"] != "DECOY_" or policy["peptide_length_range_for_occurrence_bound"] != [8, 15]:
        raise ValueError("frozen decoy prefix or peptide bound changed")
    if set(lock["normalization"]["allowed_amino_acids"]) != ALLOWED_AA:
        raise ValueError("source-lock amino-acid alphabet differs from builder")
    if lock["candidate_peptides_enumerated"] is not False:
        raise ValueError("catalog build must not enumerate peptides")
    return lock, lock_sha


def build(output_dir: Path, input_root: Path) -> dict[str, object]:
    lock, lock_sha = load_lock()
    output_dir = output_dir.resolve()
    input_root = input_root.resolve()
    if not input_root.is_dir():
        raise ValueError(f"missing input root: {input_root}")
    for ancestor in (WORKSPACE, *WORKSPACE.parents):
        if (ancestor / ".git").exists():
            if output_dir.is_relative_to(ancestor):
                raise ValueError(f"catalog output must be outside the Git checkout: {output_dir}")
            break
    if output_dir.exists():
        raise FileExistsError(f"refusing to replace existing catalog: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    components: dict[str, list[Record]] = {}
    filters: dict[str, dict[str, int]] = {}
    input_paths: dict[str, Path] = {}
    for source in SOURCE_ORDER:
        rule = lock["sources"][source]
        path = (input_root / rule["input"]).resolve()
        if not path.is_relative_to(input_root) or not path.is_file():
            raise ValueError(f"missing or out-of-root input for {source}: {path}")
        actual_sha = sha256_file(path)
        actual_bytes = path.stat().st_size
        if actual_sha != rule["sha256"] or actual_bytes != rule["bytes"]:
            raise ValueError(f"frozen input bytes differ for {source}: {path}")
        input_paths[source] = path
        components[source], filters[source] = read_component(
            path, source, rule["selection"], rule["minimum_length_aa"]
        )
        if filters[source]["raw_records"] != rule["expected_raw_records"]:
            raise ValueError(f"raw {source} count differs from source lock")
        if len(components[source]) != rule["expected_selected_records"]:
            raise ValueError(f"selected {source} count differs from source lock")

    all_records = [record for source in SOURCE_ORDER for record in components[source]]
    by_sequence: dict[str, list[Record]] = defaultdict(list)
    for record in all_records:
        by_sequence[record.sequence].append(record)
    representative_by_sequence = {
        sequence: min(
            aliases,
            key=lambda item: (SOURCE_RANK[item.source_class], item.source_accession, item.original_header),
        )
        for sequence, aliases in by_sequence.items()
    }
    representatives = sorted(
        representative_by_sequence.values(),
        key=lambda item: (SOURCE_RANK[item.source_class], item.source_accession, item.sequence_sha256),
    )
    final_accession = {
        record.sequence: f"{record.source_class}|{record.source_accession}" for record in representatives
    }
    if len(set(final_accession.values())) != len(final_accession):
        raise ValueError("duplicate merged FASTA target accessions after prefixing")

    temp_dir = Path(tempfile.mkdtemp(prefix=".phase2-catalog-build-", dir=output_dir.parent))
    try:
        components_dir = temp_dir / "components"
        components_dir.mkdir()
        combined = temp_dir / "combined-targets.fasta"
        provenance = temp_dir / "protein-provenance.tsv"
        write_fasta(
            combined,
            (
                (
                    final_accession[record.sequence],
                    f"source_class={record.source_class} sequence_sha256={record.sequence_sha256} "
                    f"source_mappings={len(by_sequence[record.sequence])}",
                    record.sequence,
                )
                for record in representatives
            ),
        )
        provenance_fields = [
            "final_accession", "final_source_class", "source_class", "source_accession",
            "original_header", "sequence_sha256", "sequence_length", "is_representative",
            "source_rank", "selection_reason",
        ]
        provenance_rows = []
        for record in sorted(
            all_records,
            key=lambda item: (
                final_accession[item.sequence], SOURCE_RANK[item.source_class],
                item.source_accession, item.original_header,
            ),
        ):
            representative = representative_by_sequence[record.sequence]
            provenance_rows.append({
                "final_accession": final_accession[record.sequence],
                "final_source_class": representative.source_class,
                "source_class": record.source_class,
                "source_accession": record.source_accession,
                "original_header": record.original_header,
                "sequence_sha256": record.sequence_sha256,
                "sequence_length": len(record.sequence),
                "is_representative": "1" if record is representative else "0",
                "source_rank": SOURCE_RANK[record.source_class],
                "selection_reason": record.selection_reason,
            })
        write_tsv(provenance, provenance_fields, provenance_rows)
        for source in SOURCE_ORDER:
            write_fasta(
                components_dir / COMPONENT_FILENAMES[source],
                (
                    (
                        f"{record.source_class}|{record.source_accession}",
                        f"source_class={record.source_class} sequence_sha256={record.sequence_sha256}",
                        record.sequence,
                    )
                    for record in sorted(
                        components[source], key=lambda item: (item.source_accession, item.original_header)
                    )
                ),
            )

        write_tsv(
            temp_dir / "input-sha256.tsv",
            ["source_class", "sha256", "bytes", "relative_path", "origin"],
            (
                {
                    "source_class": source,
                    "sha256": lock["sources"][source]["sha256"],
                    "bytes": input_paths[source].stat().st_size,
                    "relative_path": lock["sources"][source]["input"],
                    "origin": lock["sources"][source]["origin"],
                }
                for source in SOURCE_ORDER
            ),
        )
        shutil.copyfile(LOCK_PATH, temp_dir / LOCK_PATH.name)
        shutil.copyfile(Path(__file__), temp_dir / Path(__file__).name)

        sorted_sequences = sorted(record.sequence for record in representatives)
        sorted_sequence_hashes = sorted(record.sequence_sha256 for record in representatives)
        cross_source_duplicates = sum(
            len({record.source_class for record in aliases}) > 1 for aliases in by_sequence.values()
        )
        deduplication = {
            "selected_source_mappings": len(all_records),
            "unique_full_sequences": len(representatives),
            "duplicate_source_mappings_removed_from_combined_fasta": len(all_records) - len(representatives),
            "cross_source_duplicate_sequences": cross_source_duplicates,
            "combined_target_residues": sum(len(record.sequence) for record in representatives),
            "combined_raw_8_15_substring_occurrence_upper_bound": sum(
                raw_substring_bound(len(record.sequence)) for record in representatives
            ),
        }
        manifest = {
            "schema": "numivivo.melanoma.phase2.candidate_blind_search_catalog.v1",
            "status": "catalog_ready_no_raw_search",
            "candidate_peptides_enumerated": False,
            "source_lock_sha256": lock_sha,
            "builder_sha256": sha256_file(Path(__file__)),
            "normalization_reference_builder_sha256": REFERENCE_BUILDER_SHA256,
            "source_order": list(SOURCE_ORDER),
            "source_rules": lock["sources"],
            "normalization": lock["normalization"],
            "deduplication_rule": lock["deduplication"],
            "search_contract": lock["search_contract"],
            "sources": {
                source: {**source_summary(components[source]), "filter_details": filters[source]}
                for source in SOURCE_ORDER
            },
            "deduplication": deduplication,
            "content_digests": {
                "combined_fasta_sha256": sha256_file(combined),
                "protein_provenance_tsv_sha256": sha256_file(provenance),
                "sorted_unique_sequences_newline_sha256": hashlib.sha256(
                    "\n".join(sorted_sequences).encode("ascii")
                ).hexdigest(),
                "sorted_sequence_sha256_values_newline_sha256": hashlib.sha256(
                    "\n".join(sorted_sequence_hashes).encode("ascii")
                ).hexdigest(),
            },
            "validations": {
                "frozen_input_bytes_match": True,
                "frozen_source_counts_match": True,
                "sequence_alphabet_validated": True,
                "full_sequence_deduplication_applied": True,
                "every_selected_source_mapping_in_provenance": len(provenance_rows) == len(all_records),
                "no_candidate_peptide_enumeration": True,
            },
        }
        manifest_path = temp_dir / "CATALOG_MANIFEST.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        files_to_hash = [
            combined, provenance, manifest_path, temp_dir / "input-sha256.tsv",
            temp_dir / LOCK_PATH.name, temp_dir / Path(__file__).name,
        ] + [components_dir / COMPONENT_FILENAMES[source] for source in SOURCE_ORDER]
        write_tsv(
            temp_dir / "output-sha256.tsv",
            ["sha256", "bytes", "relative_path"],
            (
                {
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                    "relative_path": str(path.relative_to(temp_dir)),
                }
                for path in sorted(files_to_hash, key=lambda item: str(item.relative_to(temp_dir)))
            ),
        )
        (temp_dir / "BUILD_COMPLETE").write_text(
            f"catalog_manifest_sha256\t{sha256_file(manifest_path)}\n",
            encoding="ascii",
        )
        os.replace(temp_dir, output_dir)
    finally:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
    return deduplication


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.output_dir, args.input_root), sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"phase-2 catalog build failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
