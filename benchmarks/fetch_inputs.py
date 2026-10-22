"""Fetch the benchmark structures and their RCSB metadata.

Usage::

    python benchmarks/fetch_inputs.py [--set channels|tunnels|pockets|all] [--work-dir DIR]

Structures are downloaded with ``crevice fetch`` (deposited entry or the
listed biological assembly, mmCIF) into ``<work-dir>/pdb``; the RCSB entry,
assembly and polymer-entity metadata are saved as JSON in
``<work-dir>/rcsb`` so that the curation can be checked offline. Nothing is
written inside the repository.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INPUTS, crevice_env, default_work_dir, python_for_crevice, read_csv, run_logged  # noqa: E402

SETS = {"channels": "channels.csv", "tunnels": "tunnels.csv", "pockets": "pockets.csv"}


def structure_file(work: Path, pdb_id: str, assembly: str) -> Path:
    name = f"{pdb_id}.cif" if assembly in ("", "deposited") else f"{pdb_id}-assembly{assembly}.cif"
    return work / "pdb" / name


def rcsb_json(url: str, path: Path) -> None:
    if path.exists():
        return
    with urllib.request.urlopen(url, timeout=60) as response:
        data = json.loads(response.read().decode())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1) + "\n")


def fetch_row(work: Path, row: dict, log: Path) -> None:
    pdb_id, assembly = row["pdb_id"], row.get("assembly", "deposited")
    target = structure_file(work, pdb_id, assembly)
    if not target.exists():
        cmd = [python_for_crevice(), "-m", "crevice", "fetch", pdb_id, "--format", "cif",
               "--cache-dir", str(work / "pdb")]
        if assembly not in ("", "deposited"):
            cmd += ["--assembly", assembly]
        code, _ = run_logged(cmd, log, env=crevice_env())
        if code:
            raise SystemExit(f"fetch failed for {pdb_id}; see {log}")
    if assembly not in ("", "deposited"):
        # The assembly file names the entry XXXX; the deposited file is needed to check it.
        deposited = structure_file(work, pdb_id, "deposited")
        if not deposited.exists():
            run_logged([python_for_crevice(), "-m", "crevice", "fetch", pdb_id, "--format", "cif",
                        "--cache-dir", str(work / "pdb")], log, env=crevice_env())
    rcsb = work / "rcsb"
    rcsb_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", rcsb / f"{pdb_id}_entry.json")
    entry = json.loads((rcsb / f"{pdb_id}_entry.json").read_text())
    for asm in entry["rcsb_entry_container_identifiers"].get("assembly_ids", []):
        rcsb_json(f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{asm}", rcsb / f"{pdb_id}_assembly{asm}.json")
    for ent in entry["rcsb_entry_container_identifiers"].get("polymer_entity_ids", []):
        rcsb_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent}",
                  rcsb / f"{pdb_id}_polymer_entity{ent}.json")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--set", choices=[*SETS, "all"], default="all")
    ap.add_argument("--work-dir", type=Path, default=default_work_dir())
    ap.add_argument("--only", nargs="*", help="Entry ids to fetch (default: all in the set)")
    args = ap.parse_args(argv)
    names = list(SETS) if args.set == "all" else [args.set]
    for name in names:
        path = INPUTS / SETS[name]
        if not path.exists():
            continue
        for row in read_csv(path):
            if args.only and row["id"] not in args.only:
                continue
            fetch_row(args.work_dir, row, args.work_dir / "logs" / "fetch.log")
            print(f"{name}: {row['id']} ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
