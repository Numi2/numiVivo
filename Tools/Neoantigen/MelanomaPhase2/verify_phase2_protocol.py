#!/usr/bin/env python3
"""Read-only preflight for the frozen PXD004894 phase-2 discovery inputs.

This verifies source identity and cohort separation. It never opens search results
and cannot establish peptide presentation or a biological result.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import stat
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXPECTED = {
    "pxd_inventory": "686dbc1bcd43bf5046c57afb3627670c3b79bbc3f7b2e311ac346af2d5af4661",
    "pxd_donors": "3b258f9740fdc11fe93f7f8ef36aa8f183cdfb66a8b328f068e4bb564a1f3e96",
    "msv_inventory": "172175ffd1411a58d42bb4ddf496e5501745a61790948ccfdbd5fb91b18f4b92",
    "catalog": "a73ea6ab3e8c28262d5effec3a0baf5d98bab25153cf6da5d3b68a435db5193d",
    "provenance": "381e001d4f992a7dd42b2ed8f4e8f6db0e897f70ee47bab5921434df579d4e0c",
    "catalog_manifest": "444a4b589f1e317a5d97a07c25b64bc639200ed891c5c3cbb6d5f6f48ad8b9b9",
    "trc_bed": "e05acb0ab51094f0ca96f35ff37864f6c2addfac22e004f09712e93af8634f6f",
    "trc_fna": "558ee88208b60353312f23f6f3459f1ac3352d13de54353b8bcc86b5c34d0d7c",
    "trc_faa": "521d95a14c49fd625b33c704452430d6fdd89dca2a97a8c66a57897ee3dadbe9",
}
COMET_TEMPLATE_SHA256 = "810d1563226721f427d55261f46920ceced308b06d2b5f9b9efc5f50af6ac2fd"
SAGE_TEMPLATE_SHA256 = "73ceb13f07926c2ced63eda02f14464b87c5741f46bc3f25df552aab5a854c29"
CODONS = dict(zip(("".join(bases) for bases in itertools.product("TCAG", repeat=3)),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def file_sha(path):
    require(not path.is_symlink() and path.is_file() and stat.S_ISREG(path.stat().st_mode),
            f"not a regular, nonsymlink file: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta(path):
    records = {}
    key = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if key is not None:
                    require(key not in records, f"duplicate FASTA accession: {key}")
                    records[key] = "".join(chunks).upper()
                key = line[1:].split()[0]
                chunks = []
            elif line:
                require(key is not None, f"FASTA sequence before header: {path}")
                chunks.append(line)
    if key is not None:
        require(key not in records, f"duplicate FASTA accession: {key}")
        records[key] = "".join(chunks).upper()
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in EXPECTED:
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    paths = {name: getattr(args, name) for name in EXPECTED}
    hashes = {name: file_sha(path) for name, path in paths.items()}
    for name, actual in hashes.items():
        require(actual == EXPECTED[name], f"frozen {name} SHA-256 mismatch: {actual}")
    require(file_sha(HERE / "comet_txt_only.params.template") == COMET_TEMPLATE_SHA256,
            "frozen Comet template changed")
    sage_path = HERE / "sage_production.json.template"
    require(file_sha(sage_path) == SAGE_TEMPLATE_SHA256, "frozen Sage template changed")
    sage = json.loads(sage_path.read_text())
    require(sage["database"]["fasta"] == "__CATALOG_FASTA__" and
            sage["output_directory"] == "__OUTPUT_DIRECTORY__" and
            sage["mzml_paths"] == ["__MZML_PATH__"] and
            sage["database"]["generate_decoys"] is True and
            sage["database"]["decoy_tag"] == "rev_" and sage["write_pin"] is False,
            "Sage template must render only the frozen catalog, output, and one mzML")

    pxd = tsv(paths["pxd_inventory"])
    donors = tsv(paths["pxd_donors"])
    msv = tsv(paths["msv_inventory"])
    require(len(pxd) == 88 and len(donors) == 25 and len(msv) == 16,
            "unexpected cohort row count")
    require(all(row["dataset_accession"] == "PXD004894" and row["raw_file"].lower().endswith(".raw")
                and row["assay_name"] == row["raw_file"][:-4] and
                ("HLAp" in row["assay_name"] or "HLA-I-p" in row["assay_name"])
                for row in pxd), "PXD row is not a pinned HLA-I RAW")
    require(len({row["raw_file"] for row in pxd}) == 88 and
            len({row["pride_file_accession"] for row in pxd}) == 88, "PXD file duplication")
    by_donor = collections.defaultdict(list)
    for row in pxd:
        by_donor[row["donor"]].append(row)
    require(len(by_donor) == 25 and {row["donor"] for row in donors} == set(by_donor),
            "PXD donor summary mismatch")
    for row in donors:
        runs = by_donor[row["donor"]]
        require(int(row["hla_i_raw_count"]) == len(runs) and
                int(row["total_raw_bytes"]) == sum(int(run["raw_bytes"]) for run in runs) and
                {run["source_name"] for run in runs} == {row["source_name"]} and
                {run["hla_typing_two_field"] for run in runs} == {row["hla_typing_two_field"]} and
                row["mapping_status"] == "complete_one_donor_one_source_one_typing_all_runs",
                f"PXD donor aggregation changed: {row['donor']}")
    heldout = collections.Counter(row["donor"] for row in msv)
    require(heldout == {"Mel02": 5, "Mel11": 11} and
            all(row["dataset_accession"] == "MSV000084787" for row in msv) and
            len({row["mzml_file"] for row in msv}) == 16 and
            not set(heldout).intersection(by_donor), "held-out cohort overlap or mismatch")

    manifest = json.loads(paths["catalog_manifest"].read_text())
    require(manifest["schema"] == "numivivo.melanoma.phase2.candidate_blind_search_catalog.v1" and
            manifest["candidate_peptides_enumerated"] is False and
            manifest["source_order"] == ["CON", "UPR", "TRC"] and
            manifest["search_contract"]["decoy_search"] == 1 and
            manifest["content_digests"]["combined_fasta_sha256"] == hashes["catalog"] and
            manifest["content_digests"]["protein_provenance_tsv_sha256"] == hashes["provenance"],
            "catalog lock is inconsistent")
    source_counts = {name: manifest["sources"][name]["selected_records"] for name in manifest["source_order"]}
    require(source_counts == {"CON": 264, "UPR": 20652, "TRC": 6771}, "catalog source counts changed")
    catalog = fasta(paths["catalog"])
    provenance = tsv(paths["provenance"])
    require(len(catalog) == 27596 and len(provenance) == 27687, "catalog size changed")
    require(collections.Counter(row["source_class"] for row in provenance) == source_counts,
            "catalog provenance counts changed")
    require(all(row["final_accession"] in catalog and
                hashlib.sha256(catalog[row["final_accession"]].encode()).hexdigest() == row["sequence_sha256"]
                for row in provenance), "catalog provenance sequence mismatch")

    faa = fasta(paths["trc_faa"])
    fna = fasta(paths["trc_fna"])
    bed = {}
    with paths["trc_bed"].open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) >= 12 and fields[3] not in bed, "invalid or duplicate TRC BED12 row")
            bed[fields[3]] = fields
    require(len(faa) == 10127 and set(faa) == set(fna) == set(bed), "TRC source accession mismatch")
    trc_selected = {row["source_accession"] for row in provenance if row["source_class"] == "TRC"}
    require(len(trc_selected) == 6771 and trc_selected <= set(faa) and
            trc_selected == {key for key, sequence in faa.items() if len(sequence) >= 8},
            "TRC selection changed")
    for key, protein in faa.items():
        fields = bed[key]
        blocks = [int(value) for value in fields[10].rstrip(",").split(",")]
        require(len(fna[key]) == 3 * (len(protein) + 1) and
                len(blocks) == int(fields[9]) and sum(blocks) == len(fna[key]) and
                int(fields[2]) > int(fields[1]) and fields[5] in ("+", "-"),
                f"TRC spliced frame/BED length mismatch: {key}")
        translated = "".join(CODONS.get(fna[key][offset:offset + 3], "?")
                             for offset in range(0, len(fna[key]), 3))
        require(protein.startswith("M") and translated[1:-1] == protein[1:] and
                translated[-1] == "*" and "*" not in translated[:-1] and "?" not in translated,
                f"TRC FNA/FAA translation mismatch: {key}")
    print(json.dumps({"status": "frozen_input_preflight_pass", "pxd_raws": len(pxd),
                      "pxd_donors": len(by_donor), "heldout_runs": len(msv),
                      "heldout_donors": dict(heldout), "catalog_proteins": len(catalog),
                      "catalog_source_mappings": len(provenance), "trc_primary_orfs": len(faa),
                      "trc_selected_orfs": len(trc_selected), "sha256": hashes,
                      "comet_template_sha256": COMET_TEMPLATE_SHA256,
                      "sage_template_sha256": SAGE_TEMPLATE_SHA256}, sort_keys=True))


if __name__ == "__main__":
    main()
