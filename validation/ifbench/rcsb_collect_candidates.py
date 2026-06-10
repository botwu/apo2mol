#!/usr/bin/env python3
"""Collect independent IFBench apo/holo candidate triples from RCSB PDB.

This script builds a candidate CSV only. It does not decide benchmark quality;
all generated candidates must still pass build_manifest.py and leakage_check.py.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
if str(VALIDATION_ROOT) not in sys.path:
    sys.path.insert(0, str(VALIDATION_ROOT))

from ifbench.schema import stable_case_id


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core"

WATER_AND_IONS = {
    "HOH",
    "WAT",
    "DOD",
    "H2O",
    "NA",
    "K",
    "CL",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "CO",
    "NI",
    "CD",
    "HG",
}

COMMON_CRYSTAL_ADDITIVES = {
    "ACT",
    "ACY",
    "BME",
    "BOG",
    "CIT",
    "DMS",
    "DTT",
    "EDO",
    "FMT",
    "GOL",
    "HEP",
    "IPA",
    "MES",
    "MPD",
    "PEG",
    "PG4",
    "PO4",
    "SO4",
    "TLA",
    "TRS",
}


@dataclass
class ChainMeta:
    entry_id: str
    entity_id: str
    label_asym_id: str
    auth_asym_id: str
    uniprot_ids: list[str] = field(default_factory=list)
    cluster_key: str | None = None


@dataclass
class LigandMeta:
    entry_id: str
    entity_id: str
    comp_id: str
    label_asym_id: str | None
    auth_asym_id: str | None
    formula_weight: float | None = None


@dataclass
class EntryMeta:
    entry_id: str
    title: str | None
    release_date: str | None
    resolution: float | None
    method: str | None
    protein_chains: list[ChainMeta]
    ligands: list[LigandMeta]


def request_json(url: str, payload: dict[str, Any] | None = None, timeout: int = 60, retries: int = 2) -> dict[str, Any]:
    data = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers, method="POST" if payload is not None else "GET")
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8")
                if not text.strip():
                    return {}
                return json.loads(text)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Request failed for {url}: {last_error}")


def terminal(attribute: str, operator: str, value: Any, service: str = "text") -> dict[str, Any]:
    return {
        "type": "terminal",
        "service": service,
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


def group(operator: str, nodes: list[dict[str, Any]]) -> dict[str, Any]:
    return {"type": "group", "logical_operator": operator, "nodes": nodes}


def rcsb_entry_query(kind: str, max_resolution: float, release_after: str | None) -> dict[str, Any]:
    nodes = [
        terminal("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
        terminal("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution),
        terminal("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
    ]
    if release_after:
        nodes.append(terminal("rcsb_accession_info.initial_release_date", "greater_or_equal", release_after))
    if kind == "holo":
        nodes.append(terminal("rcsb_entry_info.nonpolymer_entity_count", "greater_or_equal", 1))
    elif kind == "apo":
        nodes.append(terminal("rcsb_entry_info.nonpolymer_entity_count", "equals", 0))
    else:
        raise ValueError(f"Unknown query kind: {kind}")
    return group("and", nodes)


def search_entries(kind: str, max_entries: int, max_resolution: float, release_after: str | None, page_size: int = 1000) -> list[str]:
    entry_ids: list[str] = []
    start = 0
    while len(entry_ids) < max_entries:
        payload = {
            "query": rcsb_entry_query(kind, max_resolution=max_resolution, release_after=release_after),
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": min(page_size, max_entries - len(entry_ids))},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "desc"}],
            },
        }
        data = request_json(SEARCH_URL, payload=payload)
        rows = data.get("result_set", [])
        if not rows:
            break
        entry_ids.extend(str(row["identifier"]).lower() for row in rows if row.get("identifier"))
        start += len(rows)
        if len(rows) < page_size:
            break
    return entry_ids[:max_entries]


def uniprot_query(accession: str, apo_only: bool, max_resolution: float) -> dict[str, Any]:
    nodes = [
        terminal(
            "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "exact_match",
            accession,
        ),
        terminal("rcsb_entry_info.polymer_entity_count_protein", "greater_or_equal", 1),
        terminal("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution),
        terminal("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
    ]
    if apo_only:
        nodes.append(terminal("rcsb_entry_info.nonpolymer_entity_count", "equals", 0))
    return group("and", nodes)


def search_entries_by_uniprot(
    accessions: list[str],
    max_entries_per_uniprot: int,
    max_total_entries: int,
    max_resolution: float,
    apo_only: bool = True,
) -> list[str]:
    found: list[str] = []
    seen: set[str] = set()
    for accession in accessions:
        if len(found) >= max_total_entries:
            break
        payload = {
            "query": uniprot_query(accession, apo_only=apo_only, max_resolution=max_resolution),
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": max_entries_per_uniprot},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "asc"}],
            },
        }
        try:
            data = request_json(SEARCH_URL, payload=payload)
        except RuntimeError as exc:
            print(f"[WARN] UniProt search failed for {accession}: {exc}", file=sys.stderr)
            continue
        rows = data.get("result_set", [])
        print(f"[INFO] UniProt {accession}: {len(rows)} apo candidate entries", file=sys.stderr)
        for row in rows:
            entry_id = str(row.get("identifier", "")).lower()
            if entry_id and entry_id not in seen:
                seen.add(entry_id)
                found.append(entry_id)
                if len(found) >= max_total_entries:
                    break
    return found


def search_payload(kind: str, max_entries: int, max_resolution: float, release_after: str | None) -> dict[str, Any]:
    return {
        "query": rcsb_entry_query(kind, max_resolution=max_resolution, release_after=release_after),
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_entries},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "desc"}],
        },
    }


def read_entry_list(path: Path | None) -> list[str]:
    if not path:
        return []
    return [line.strip().lower() for line in path.read_text().splitlines() if line.strip() and not line.startswith("#")]


def first_resolution(entry: dict[str, Any]) -> float | None:
    values = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    clean = [float(value) for value in values if value is not None]
    return min(clean) if clean else None


def first_method(entry: dict[str, Any]) -> str | None:
    methods = [item.get("method") for item in entry.get("exptl", []) if item.get("method")]
    return ";".join(methods) if methods else None


def first_release_date(entry: dict[str, Any]) -> str | None:
    value = entry.get("rcsb_accession_info", {}).get("initial_release_date")
    if not value:
        return None
    return str(value).split("T", 1)[0]


def sequence_hash(value: str | None) -> str | None:
    if not value:
        return None
    return hashlib.sha1(value.encode("utf-8")).hexdigest()[:12]


def fetch_entry(entry_id: str) -> dict[str, Any]:
    return request_json(f"{DATA_URL}/entry/{entry_id}")


def fetch_polymer_entity(entry_id: str, entity_id: str) -> dict[str, Any]:
    return request_json(f"{DATA_URL}/polymer_entity/{entry_id}/{entity_id}")


def fetch_nonpolymer_entity(entry_id: str, entity_id: str) -> dict[str, Any]:
    return request_json(f"{DATA_URL}/nonpolymer_entity/{entry_id}/{entity_id}")


def chain_meta_from_dict(item: dict[str, Any]) -> ChainMeta:
    return ChainMeta(
        entry_id=str(item["entry_id"]),
        entity_id=str(item["entity_id"]),
        label_asym_id=str(item["label_asym_id"]),
        auth_asym_id=str(item["auth_asym_id"]),
        uniprot_ids=[str(value) for value in item.get("uniprot_ids", [])],
        cluster_key=item.get("cluster_key"),
    )


def ligand_meta_from_dict(item: dict[str, Any]) -> LigandMeta:
    formula_weight = item.get("formula_weight")
    if formula_weight not in (None, ""):
        formula_weight = float(formula_weight)
    return LigandMeta(
        entry_id=str(item["entry_id"]),
        entity_id=str(item["entity_id"]),
        comp_id=str(item["comp_id"]),
        label_asym_id=item.get("label_asym_id"),
        auth_asym_id=item.get("auth_asym_id"),
        formula_weight=formula_weight,
    )


def entry_meta_from_dict(item: dict[str, Any]) -> EntryMeta:
    return EntryMeta(
        entry_id=str(item["entry_id"]),
        title=item.get("title"),
        release_date=item.get("release_date"),
        resolution=float(item["resolution"]) if item.get("resolution") not in (None, "") else None,
        method=item.get("method"),
        protein_chains=[chain_meta_from_dict(value) for value in item.get("protein_chains", [])],
        ligands=[ligand_meta_from_dict(value) for value in item.get("ligands", [])],
    )


def load_metadata_cache(path: Path | None) -> dict[str, EntryMeta]:
    if not path or not path.exists():
        return {}
    payload = json.loads(path.read_text())
    records: dict[str, EntryMeta] = {}
    if "entries" in payload:
        values = payload["entries"].values()
    else:
        values = list(payload.get("holo", [])) + list(payload.get("apo", []))
    for item in values:
        record = entry_meta_from_dict(item)
        records[record.entry_id] = record
    return records


def dump_metadata_cache(path: Path, records: dict[str, EntryMeta]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {"entries": {entry_id: record.__dict__ for entry_id, record in sorted(records.items())}},
            default=lambda obj: obj.__dict__,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def filter_ligands_by_weight(record: EntryMeta, min_ligand_formula_weight: float) -> EntryMeta:
    record.ligands = [
        ligand
        for ligand in record.ligands
        if ligand.formula_weight is None or ligand.formula_weight >= min_ligand_formula_weight
    ]
    return record


def normalize_formula_weight(value: Any) -> float | None:
    if value is None or value == "":
        return None
    weight = float(value)
    # RCSB nonpolymer formula_weight is reported in kDa for many entries.
    if 0.0 < weight < 5.0:
        weight *= 1000.0
    return weight


def entity_ids(entry: dict[str, Any], key: str) -> list[str]:
    values = entry.get("rcsb_entry_container_identifiers", {}).get(key) or []
    return [str(value) for value in values]


def parse_polymer_chains(entry_id: str, entry: dict[str, Any]) -> list[ChainMeta]:
    chains: list[ChainMeta] = []
    for entity_id in entity_ids(entry, "polymer_entity_ids"):
        try:
            entity = fetch_polymer_entity(entry_id, entity_id)
        except RuntimeError as exc:
            print(f"[WARN] skip polymer entity {entry_id}/{entity_id}: {exc}", file=sys.stderr)
            continue
        entity_poly = entity.get("entity_poly", {})
        if str(entity_poly.get("rcsb_entity_polymer_type", "")).lower() != "protein":
            continue
        ids = entity.get("rcsb_polymer_entity_container_identifiers", {})
        asym_ids = [str(value) for value in ids.get("asym_ids") or []]
        auth_asym_ids = [str(value) for value in ids.get("auth_asym_ids") or []]
        uniprot_ids = [str(value) for value in ids.get("uniprot_ids") or []]
        cluster_key = uniprot_ids[0] if uniprot_ids else sequence_hash(entity_poly.get("pdbx_seq_one_letter_code_can"))
        for idx, label_asym in enumerate(asym_ids):
            auth_asym = auth_asym_ids[idx] if idx < len(auth_asym_ids) else label_asym
            chains.append(
                ChainMeta(
                    entry_id=entry_id,
                    entity_id=entity_id,
                    label_asym_id=label_asym,
                    auth_asym_id=auth_asym,
                    uniprot_ids=uniprot_ids,
                    cluster_key=cluster_key,
                )
            )
    return chains


def parse_ligands(entry_id: str, entry: dict[str, Any], include_additives: bool, min_ligand_formula_weight: float) -> list[LigandMeta]:
    ligands: list[LigandMeta] = []
    excluded = set(WATER_AND_IONS)
    if not include_additives:
        excluded.update(COMMON_CRYSTAL_ADDITIVES)
    for entity_id in entity_ids(entry, "non_polymer_entity_ids"):
        try:
            entity = fetch_nonpolymer_entity(entry_id, entity_id)
        except RuntimeError as exc:
            print(f"[WARN] skip nonpolymer entity {entry_id}/{entity_id}: {exc}", file=sys.stderr)
            continue
        ids = entity.get("rcsb_nonpolymer_entity_container_identifiers", {})
        comp = (entity.get("nonpolymer_comp", {}) or {}).get("chem_comp", {}) or {}
        rcsb_nonpoly = entity.get("rcsb_nonpolymer_entity", {}) or {}
        pdbx_nonpoly = entity.get("pdbx_entity_nonpoly", {}) or {}
        comp_id = str(
            comp.get("id")
            or pdbx_nonpoly.get("comp_id")
            or ids.get("nonpolymer_comp_id")
            or ""
        ).upper()
        if not comp_id or comp_id in excluded:
            continue
        asym_ids = [str(value) for value in ids.get("asym_ids") or []]
        auth_asym_ids = [str(value) for value in ids.get("auth_asym_ids") or []]
        formula_weight = normalize_formula_weight(comp.get("formula_weight") or rcsb_nonpoly.get("formula_weight"))
        if formula_weight is not None and formula_weight < min_ligand_formula_weight:
            continue
        for idx, label_asym in enumerate(asym_ids or [None]):
            auth_asym = auth_asym_ids[idx] if idx < len(auth_asym_ids) else label_asym
            ligands.append(
                LigandMeta(
                    entry_id=entry_id,
                    entity_id=entity_id,
                    comp_id=comp_id,
                    label_asym_id=label_asym,
                    auth_asym_id=auth_asym,
                    formula_weight=formula_weight,
                )
            )
    return ligands


def fetch_entry_meta(entry_id: str, include_additives: bool, min_ligand_formula_weight: float) -> EntryMeta | None:
    try:
        entry = fetch_entry(entry_id)
    except RuntimeError as exc:
        print(f"[WARN] skip entry {entry_id}: {exc}", file=sys.stderr)
        return None
    return EntryMeta(
        entry_id=entry_id,
        title=(entry.get("struct", {}) or {}).get("title"),
        release_date=first_release_date(entry),
        resolution=first_resolution(entry),
        method=first_method(entry),
        protein_chains=parse_polymer_chains(entry_id, entry),
        ligands=parse_ligands(
            entry_id,
            entry,
            include_additives=include_additives,
            min_ligand_formula_weight=min_ligand_formula_weight,
        ),
    )


def collect_metadata(
    entry_ids: list[str],
    include_additives: bool,
    min_ligand_formula_weight: float,
    cache: dict[str, EntryMeta] | None = None,
    cache_path: Path | None = None,
) -> list[EntryMeta]:
    records: list[EntryMeta] = []
    cache = cache if cache is not None else {}
    for idx, entry_id in enumerate(entry_ids, start=1):
        if entry_id in cache:
            print(f"[{idx}/{len(entry_ids)}] metadata {entry_id} (cache)", file=sys.stderr)
            records.append(filter_ligands_by_weight(cache[entry_id], min_ligand_formula_weight))
            continue
        print(f"[{idx}/{len(entry_ids)}] metadata {entry_id}", file=sys.stderr)
        record = fetch_entry_meta(
            entry_id,
            include_additives=include_additives,
            min_ligand_formula_weight=min_ligand_formula_weight,
        )
        if record is not None and record.protein_chains:
            records.append(record)
            cache[record.entry_id] = record
            if cache_path:
                dump_metadata_cache(cache_path, cache)
    return records


def row_paths(holo: EntryMeta, apo: EntryMeta, ligand: LigandMeta, raw_prefix: str) -> tuple[str, str, str]:
    ligand_asym = ligand.label_asym_id or ligand.auth_asym_id or ligand.entity_id
    ligand_name = f"{holo.entry_id}_{ligand.comp_id}_{ligand_asym}.sdf"
    return (
        f"{raw_prefix}/structures/{holo.entry_id}.cif",
        f"{raw_prefix}/structures/{apo.entry_id}.cif",
        f"{raw_prefix}/ligands/{ligand_name}",
    )


def write_metadata_cache(path: Path, holo: list[EntryMeta], apo: list[EntryMeta]) -> None:
    records = {record.entry_id: record for record in holo + apo}
    dump_metadata_cache(path, records)


def make_candidate_rows(
    holo_records: list[EntryMeta],
    apo_records: list[EntryMeta],
    raw_prefix: str,
    max_pairs_per_holo: int,
) -> list[dict[str, Any]]:
    apo_by_key: dict[str, list[tuple[EntryMeta, ChainMeta]]] = {}
    for apo in apo_records:
        if apo.ligands:
            continue
        for chain in apo.protein_chains:
            if chain.cluster_key:
                apo_by_key.setdefault(chain.cluster_key, []).append((apo, chain))

    rows: list[dict[str, Any]] = []
    for holo in holo_records:
        for holo_chain in holo.protein_chains:
            if not holo_chain.cluster_key:
                continue
            apo_options = [item for item in apo_by_key.get(holo_chain.cluster_key, []) if item[0].entry_id != holo.entry_id]
            if not apo_options:
                continue
            for ligand in holo.ligands:
                for rank, (apo, apo_chain) in enumerate(apo_options[:max_pairs_per_holo], start=1):
                    holo_file, apo_file, ligand_file = row_paths(holo, apo, ligand, raw_prefix=raw_prefix)
                    case_id = stable_case_id(
                        [
                            "rcsb",
                            holo.entry_id,
                            holo_chain.label_asym_id,
                            apo.entry_id,
                            apo_chain.label_asym_id,
                            ligand.comp_id,
                            ligand.label_asym_id or "",
                        ]
                    )
                    uniprot = (holo_chain.uniprot_ids or apo_chain.uniprot_ids or [None])[0]
                    rows.append(
                        {
                            "case_id": case_id,
                            "source": "rcsb",
                            "holo_pdb_id": holo.entry_id,
                            "apo_pdb_id": apo.entry_id,
                            "holo_chain_id": holo_chain.label_asym_id,
                            "apo_chain_id": apo_chain.label_asym_id,
                            "ligand_id": ligand.comp_id,
                            "holo_structure_file": holo_file,
                            "apo_structure_file": apo_file,
                            "ligand_file": ligand_file,
                            "uniprot_id": uniprot or "",
                            "cluster_key": holo_chain.cluster_key or "",
                            "holo_release_date": holo.release_date or "",
                            "apo_release_date": apo.release_date or "",
                            "holo_resolution": holo.resolution if holo.resolution is not None else "",
                            "apo_resolution": apo.resolution if apo.resolution is not None else "",
                            "holo_method": holo.method or "",
                            "apo_method": apo.method or "",
                            "holo_entity_id": holo_chain.entity_id,
                            "apo_entity_id": apo_chain.entity_id,
                            "holo_auth_chain_id": holo_chain.auth_asym_id,
                            "apo_auth_chain_id": apo_chain.auth_asym_id,
                            "ligand_entity_id": ligand.entity_id,
                            "ligand_label_asym_id": ligand.label_asym_id or "",
                            "ligand_auth_asym_id": ligand.auth_asym_id or "",
                            "ligand_formula_weight": ligand.formula_weight if ligand.formula_weight is not None else "",
                            "candidate_rank": rank,
                            "holo_title": holo.title or "",
                            "apo_title": apo.title or "",
                        }
                    )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    if not fieldnames:
        fieldnames = [
            "case_id",
            "source",
            "holo_pdb_id",
            "apo_pdb_id",
            "holo_chain_id",
            "apo_chain_id",
            "ligand_id",
            "holo_structure_file",
            "apo_structure_file",
            "ligand_file",
        ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-candidates", required=True, type=Path)
    parser.add_argument("--holo-entry-list", type=Path, help="Optional newline-separated holo PDB IDs.")
    parser.add_argument("--apo-entry-list", type=Path, help="Optional newline-separated apo PDB IDs.")
    parser.add_argument("--max-holo-entries", type=int, default=200)
    parser.add_argument("--max-apo-entries", type=int, default=2000)
    parser.add_argument("--max-apo-per-uniprot", type=int, default=50)
    parser.add_argument("--max-pairs-per-holo", type=int, default=3)
    parser.add_argument("--max-resolution", type=float, default=2.5)
    parser.add_argument("--holo-release-after", default="2024-05-01", help="ISO date filter for independent post-Apo2Mol holo seeds; empty disables.")
    parser.add_argument(
        "--apo-search-mode",
        choices=["matched-uniprot", "broad"],
        default="matched-uniprot",
        help="matched-uniprot searches apo entries for holo UniProt IDs; broad scans recent apo entries.",
    )
    parser.add_argument("--raw-prefix", default="validation/ifbench/out/raw/rcsb")
    parser.add_argument("--metadata-cache", type=Path)
    parser.add_argument("--query-json", type=Path, help="Optional path to write the RCSB Search API payloads.")
    parser.add_argument("--dry-run", action="store_true", help="Write query JSON and exit without contacting RCSB.")
    parser.add_argument("--include-common-additives", action="store_true")
    parser.add_argument("--min-ligand-formula-weight", type=float, default=120.0)
    args = parser.parse_args()

    release_after = args.holo_release_after or None
    query_payloads = {
        "holo": search_payload("holo", args.max_holo_entries, max_resolution=args.max_resolution, release_after=release_after),
        "apo": search_payload("apo", args.max_apo_entries, max_resolution=args.max_resolution, release_after=None),
    }
    if args.query_json:
        args.query_json.parent.mkdir(parents=True, exist_ok=True)
        args.query_json.write_text(json.dumps(query_payloads, indent=2, sort_keys=True) + "\n")
    if args.dry_run:
        write_csv(args.out_candidates, [])
        print(
            json.dumps(
                {
                    "dry_run": True,
                    "out_candidates": str(args.out_candidates),
                    "query_json": str(args.query_json) if args.query_json else None,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return

    holo_ids = read_entry_list(args.holo_entry_list)
    if not holo_ids:
        holo_ids = search_entries("holo", args.max_holo_entries, max_resolution=args.max_resolution, release_after=release_after)
    metadata_cache = load_metadata_cache(args.metadata_cache)
    if metadata_cache:
        print(f"Loaded {len(metadata_cache)} cached RCSB metadata records", file=sys.stderr)
    print(f"Collecting metadata for {len(holo_ids)} holo entries", file=sys.stderr)
    holo_records = collect_metadata(
        holo_ids,
        include_additives=args.include_common_additives,
        min_ligand_formula_weight=args.min_ligand_formula_weight,
        cache=metadata_cache,
        cache_path=args.metadata_cache,
    )
    apo_ids = read_entry_list(args.apo_entry_list)
    if not apo_ids:
        if args.apo_search_mode == "matched-uniprot":
            accessions = sorted(
                {
                    uniprot
                    for record in holo_records
                    for chain in record.protein_chains
                    for uniprot in chain.uniprot_ids
                    if uniprot
                }
            )
            print(f"Searching apo entries for {len(accessions)} holo UniProt IDs", file=sys.stderr)
            apo_ids = search_entries_by_uniprot(
                accessions,
                max_entries_per_uniprot=args.max_apo_per_uniprot,
                max_total_entries=args.max_apo_entries,
                max_resolution=args.max_resolution,
                apo_only=True,
            )
        else:
            apo_ids = search_entries("apo", args.max_apo_entries, max_resolution=args.max_resolution, release_after=None)
    print(f"Collecting metadata for {len(apo_ids)} apo entries", file=sys.stderr)
    apo_records = collect_metadata(
        apo_ids,
        include_additives=args.include_common_additives,
        min_ligand_formula_weight=args.min_ligand_formula_weight,
        cache=metadata_cache,
        cache_path=args.metadata_cache,
    )
    if args.metadata_cache:
        write_metadata_cache(args.metadata_cache, holo_records, apo_records)

    rows = make_candidate_rows(
        holo_records=holo_records,
        apo_records=apo_records,
        raw_prefix=args.raw_prefix.rstrip("/"),
        max_pairs_per_holo=args.max_pairs_per_holo,
    )
    write_csv(args.out_candidates, rows)
    print(
        json.dumps(
            {
                "holo_entries": len(holo_records),
                "apo_entries": len(apo_records),
                "candidate_rows": len(rows),
                "out_candidates": str(args.out_candidates),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
