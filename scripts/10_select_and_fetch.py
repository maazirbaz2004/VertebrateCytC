"""Phase 1: select and fetch one verified cytochrome c sequence per species."""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from Bio import Align, Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve  # noqa: E402

FIELDS = (
    "accession,id,gene_names,protein_name,organism_name,"
    "organism_id,length,reviewed,sequence,fragment"
)
BAD_NAME = re.compile(
    r"oxidase|c1\b|peroxidase|reductase|-type|like|testis|subunit|"
    r"chaperone|assembly|hinge|complex|synthesis",
    re.I,
)
HUMAN_ANCHOR_ACC = "P99999"


def make_headers(cfg: dict[str, Any]) -> dict[str, str]:
    """Build polite HTTP headers with contact email from config."""
    email = cfg["data_collection"]["ncbi_email"]
    return {"User-Agent": f"cytc-phylo/1.0 ({email})"}


def request_json(
    url: str,
    cfg: dict[str, Any],
    params: dict[str, Any] | None = None,
    retries: int = 3,
) -> dict[str, Any]:
    """GET JSON from UniProt with rate-limit pause and retry on transient errors."""
    pause = float(cfg["data_collection"]["request_pause_sec"])
    headers = make_headers(cfg)
    for attempt in range(1, retries + 1):
        resp = requests.get(url, params=params, headers=headers, timeout=60)
        time.sleep(pause)
        if resp.status_code in (429, 500, 502, 503, 504) and attempt < retries:
            time.sleep(2.0 * attempt)
            continue
        resp.raise_for_status()
        return resp.json()
    raise RuntimeError(f"Failed request after {retries} attempts: {url}")


def fetch_uniprot_entry(accession: str, cfg: dict[str, Any]) -> dict[str, Any]:
    """Fetch a single full UniProtKB JSON record by accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    return request_json(url, cfg)


def search_uniprot(query: str, cfg: dict[str, Any], size: int = 25) -> list[dict[str, Any]]:
    """Run a UniProtKB REST search and return candidate records."""
    url = cfg["data_collection"]["uniprot_rest"]
    params = {"query": query, "format": "json", "fields": FIELDS, "size": size}
    return request_json(url, cfg, params=params).get("results", [])


def extract_protein_name(entry: dict[str, Any]) -> str:
    """Extract recommended or submitted full protein name from a UniProt entry."""
    pd_info = entry.get("proteinDescription", {})
    rec = pd_info.get("recommendedName", {}).get("fullName", {}).get("value", "")
    if rec:
        return rec
    subs = pd_info.get("submissionNames", [])
    if subs:
        return subs[0].get("fullName", {}).get("value", "")
    return ""


def extract_gene_names(entry: dict[str, Any]) -> tuple[list[str], list[str]]:
    """Return (primary_genes, all_genes_and_synonyms) from a UniProt entry."""
    primary: list[str] = []
    all_names: list[str] = []
    for g in entry.get("genes", []):
        name = g.get("geneName", {}).get("value")
        if name:
            primary.append(name)
            all_names.append(name)
        for syn in g.get("synonyms", []):
            sval = syn.get("value")
            if sval:
                all_names.append(sval)
    return primary, all_names


def is_reviewed_entry(entry: dict[str, Any]) -> bool:
    """Return True only if UniProt entry is Swiss-Prot reviewed."""
    entry_type = entry.get("entryType", "")
    return entry_type.startswith("UniProtKB reviewed")


def validate_entry(entry: dict[str, Any], sp: pd.Series, cfg: dict[str, Any]) -> list[str]:
    """Validate a candidate record against organism, length, motif, name, and alphabet."""
    dc = cfg["data_collection"]
    std_alpha = set(cfg["sequence_prep"]["allowed_alphabet"])
    why: list[str] = []

    seq = entry.get("sequence", {}).get("value", "")
    pname = extract_protein_name(entry)
    flag = (entry.get("proteinDescription", {}).get("flag", "") or "").lower()
    is_outgroup = str(sp["role"]).strip().lower() == "outgroup"
    lo = int(dc["outgroup_len_min"] if is_outgroup else dc["vertebrate_len_min"])
    hi = int(dc["outgroup_len_max"] if is_outgroup else dc["vertebrate_len_max"])

    if "fragment" in flag:
        why.append("fragment flag present")
    if not (lo <= len(seq) <= hi):
        why.append(f"length {len(seq)} outside [{lo},{hi}]")
    if not re.search(dc["heme_motif_regex"], seq):
        why.append("no CxxCH motif")
    bad_res = sorted(set(seq) - std_alpha)
    if bad_res:
        why.append(f"non-standard residues: {''.join(bad_res)}")
    if not pname or "cytochrome c" not in pname.lower():
        why.append(f"name lacks 'cytochrome c': {pname!r}")
    elif BAD_NAME.search(pname):
        why.append(f"bad protein name: {pname!r}")

    tax = entry.get("organism", {}).get("taxonId")
    expected_tax = int(sp["taxon_id"])
    valid_yeast = sp["species"] == "Saccharomyces_cerevisiae" and tax in (559292, 4932)
    if tax is None or (int(tax) != expected_tax and not valid_yeast):
        why.append(f"taxon mismatch {tax} != {expected_tax}")
    return why


def aligned_identity(seq_a: str, seq_b: str) -> float:
    """Compute global sequence identity fraction using Bio.Align.PairwiseAligner."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    score = aligner.score(seq_a, seq_b)
    return float(score) / float(max(len(seq_a), len(seq_b), 1))


def rank_key(
    entry: dict[str, Any],
    sp: pd.Series,
    human_seq: str,
) -> tuple[Any, ...]:
    """Compute deterministic sorting key (lowest wins) for passing candidates."""
    acc = entry["primaryAccession"]
    reviewed_rank = 0 if is_reviewed_entry(entry) else 1
    hint = str(sp.get("uniprot_gene_hint") or "").strip().lower()
    primary_genes, all_genes = extract_gene_names(entry)
    genes_lower = {g.lower() for g in all_genes}

    exact_targets = {t for t in (hint, "cycs", "cycsa", "cyc-a") if t}
    if genes_lower & exact_targets:
        gene_rank = 0
    elif any(g.startswith("cyc") and g != "cyct" for g in genes_lower):
        gene_rank = 1
    else:
        gene_rank = 2

    pname = extract_protein_name(entry).lower()
    name_rank = 0 if "domain-containing" not in pname else 1

    seq = entry.get("sequence", {}).get("value", "")
    is_outgroup = str(sp["role"]).strip().lower() == "outgroup"
    len_rank = 0 if (is_outgroup or len(seq) in (104, 105)) else abs(len(seq) - 105)
    ident = aligned_identity(human_seq, seq) if human_seq else 0.0
    return (reviewed_rank, gene_rank, name_rank, len_rank, -round(ident, 6), acc)


def ncbi_fallback(sp: pd.Series, cfg: dict[str, Any], log: Any) -> list[dict[str, Any]]:
    """Fallback query to NCBI Protein via Entrez, mapped into UniProt-like dicts."""
    dc = cfg["data_collection"]
    Entrez.email = dc["ncbi_email"]
    taxid = int(sp["taxon_id"])
    term = f"cytochrome c[Title] AND txid{taxid}[Organism:exp]"
    log.info("Tier 4 NCBI fallback search for %s: %s", sp["species"], term)
    with Entrez.esearch(db="protein", term=term, retmax=15) as handle:
        ids = Entrez.read(handle).get("IdList", [])
    time.sleep(float(dc["request_pause_sec"]))
    if not ids:
        return []

    records: list[dict[str, Any]] = []
    with Entrez.efetch(db="protein", id=",".join(ids), rettype="gb", retmode="text") as handle:
        for rec in SeqIO.parse(handle, "genbank"):
            tax_val = None
            for feat in rec.features:
                if feat.type == "source":
                    for xref in feat.qualifiers.get("db_xref", []):
                        if xref.startswith("taxon:"):
                            tax_val = int(xref.split(":", 1)[1])
            records.append(
                {
                    "primaryAccession": rec.id,
                    "entryType": "NCBI GenPept",
                    "organism": {
                        "scientificName": rec.annotations.get("organism", ""),
                        "taxonId": tax_val or taxid,
                    },
                    "proteinDescription": {
                        "recommendedName": {"fullName": {"value": rec.description}}
                    },
                    "genes": [],
                    "sequence": {"value": str(rec.seq), "length": len(rec.seq)},
                    "_source_db": "NCBI",
                    "_source_url": f"https://www.ncbi.nlm.nih.gov/protein/{rec.id}",
                }
            )
    time.sleep(float(dc["request_pause_sec"]))
    return records


def collect_candidates_for_species(
    sp: pd.Series,
    cfg: dict[str, Any],
    log: Any,
) -> tuple[list[tuple[dict[str, Any], str, str]], list[dict[str, Any]]]:
    """Evaluate Tiers 1-4 in order, stopping at first tier with >= 1 passing hit."""
    dc = cfg["data_collection"]
    species_name = str(sp["species"])
    taxid = int(sp["taxon_id"])
    gene = str(sp.get("uniprot_gene_hint") or "").strip()
    cand_acc = str(sp.get("candidate_accession") or "").strip()
    if cand_acc.lower() == "nan":
        cand_acc = ""

    is_outgroup = str(sp["role"]).strip().lower() == "outgroup"
    lo = int(dc["outgroup_len_min"] if is_outgroup else dc["vertebrate_len_min"])
    hi = int(dc["outgroup_len_max"] if is_outgroup else dc["vertebrate_len_max"])

    tiers: list[tuple[str, str, list[dict[str, Any]]]] = []
    if cand_acc:
        url = f"https://rest.uniprot.org/uniprotkb/{cand_acc}.json"
        tiers.append(("Tier1_accession", url, [fetch_uniprot_entry(cand_acc, cfg)]))

    if gene:
        q_rev = f"(gene_exact:{gene}) AND (organism_id:{taxid}) AND (reviewed:true)"
        q_all = f"(gene_exact:{gene}) AND (organism_id:{taxid})"
        tiers.append(("Tier2_gene_reviewed", dc["uniprot_rest"], search_uniprot(q_rev, cfg)))
        tiers.append(("Tier2_gene_all", dc["uniprot_rest"], search_uniprot(q_all, cfg)))

    q3_rev = (
        f'(organism_id:{taxid}) AND (length:[{lo} TO {hi}]) '
        f'AND (protein_name:"cytochrome c" OR gene:{gene}) AND (reviewed:true)'
    )
    q3_all = (
        f'(organism_id:{taxid}) AND (length:[{lo} TO {hi}]) '
        f'AND (protein_name:"cytochrome c" OR gene:{gene})'
    )
    tiers.append(("Tier3_name_reviewed", dc["uniprot_rest"], search_uniprot(q3_rev, cfg)))
    tiers.append(("Tier3_name_all", dc["uniprot_rest"], search_uniprot(q3_all, cfg)))

    passing: list[tuple[dict[str, Any], str, str]] = []
    audit_rows: list[dict[str, Any]] = []
    seen_accs: set[str] = set()

    for tier_name, tier_url, hits in tiers:
        for entry in hits:
            acc = entry["primaryAccession"]
            if acc in seen_accs:
                continue
            seen_accs.add(acc)
            reasons = validate_entry(entry, sp, cfg)
            passed = len(reasons) == 0
            audit_rows.append(
                {
                    "species": species_name,
                    "tier": tier_name,
                    "accession": acc,
                    "reviewed": is_reviewed_entry(entry),
                    "organism": entry.get("organism", {}).get("scientificName", ""),
                    "taxon_id": entry.get("organism", {}).get("taxonId"),
                    "gene_names": ";".join(extract_gene_names(entry)[0]),
                    "protein_name": extract_protein_name(entry),
                    "length": len(entry.get("sequence", {}).get("value", "")),
                    "passed": passed,
                    "failure_reasons": "; ".join(reasons),
                }
            )
            if passed:
                passing.append((entry, tier_name, tier_url))
        if passing:
            log.info("%s: %d passing hit(s) at %s", species_name, len(passing), tier_name)
            return passing, audit_rows

    # Tier 4: NCBI Entrez fallback
    for entry in ncbi_fallback(sp, cfg, log):
        acc = entry["primaryAccession"]
        reasons = validate_entry(entry, sp, cfg)
        passed = len(reasons) == 0
        audit_rows.append(
            {
                "species": species_name,
                "tier": "Tier4_NCBI",
                "accession": acc,
                "reviewed": False,
                "organism": entry.get("organism", {}).get("scientificName", ""),
                "taxon_id": entry.get("organism", {}).get("taxonId"),
                "gene_names": "",
                "protein_name": extract_protein_name(entry),
                "length": len(entry.get("sequence", {}).get("value", "")),
                "passed": passed,
                "failure_reasons": "; ".join(reasons),
            }
        )
        if passed:
            passing.append((entry, "Tier4_NCBI", entry["_source_url"]))

    if not passing:
        raise RuntimeError(f"No passing candidate found for {species_name} (taxon {taxid})")
    return passing, audit_rows


def process_species(
    species_df: pd.DataFrame,
    cfg: dict[str, Any],
    human_seq: str,
    uniprot_dir: Path,
    log: Any,
) -> tuple[list[str], list[dict[str, Any]], list[dict[str, Any]], pd.DataFrame]:
    """Retrieve, rank, save raw JSON, and assemble FASTA + provenance rows."""
    fasta_lines: list[str] = []
    provenance_rows: list[dict[str, Any]] = []
    all_audit: list[dict[str, Any]] = []
    updated_df = species_df.copy()
    today_iso = date.today().isoformat()

    for idx, sp in updated_df.iterrows():
        species_name = str(sp["species"])
        passing, audit_rows = collect_candidates_for_species(sp, cfg, log)
        all_audit.extend(audit_rows)

        passing.sort(key=lambda item: rank_key(item[0], sp, human_seq))
        chosen, tier_name, _tier_url = passing[0]
        acc = chosen["primaryAccession"]
        db = chosen.get("_source_db", "UniProt")

        if db == "UniProt":
            full_json = fetch_uniprot_entry(acc, cfg)
            raw_json_path = uniprot_dir / f"{acc}.json"
            raw_json_path.write_text(json.dumps(full_json, indent=2), encoding="utf-8")
            url_used = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
            chosen = full_json
        else:
            url_used = chosen["_source_url"]

        seq = chosen["sequence"]["value"]
        primary_genes, _ = extract_gene_names(chosen)
        pname = extract_protein_name(chosen)
        reviewed = is_reviewed_entry(chosen)
        ident = aligned_identity(human_seq, seq)

        log.info(
            "CHOSEN %s -> %s (%s, reviewed=%s, len=%d, genes=%s, id_human=%.3f, from %d passing)",
            species_name,
            acc,
            db,
            reviewed,
            len(seq),
            primary_genes,
            ident,
            len(passing),
        )
        updated_df.at[idx, "candidate_accession"] = acc
        fasta_lines.append(f">{species_name}|{acc}|{db}\n{seq}\n")
        provenance_rows.append(
            {
                "species": species_name,
                "taxon_id": chosen.get("organism", {}).get("taxonId"),
                "accession": acc,
                "database": db,
                "reviewed": reviewed,
                "gene_name": ";".join(primary_genes),
                "protein_name": pname,
                "length": len(seq),
                "download_date": today_iso,
                "url": url_used,
                "notes": f"selected via {tier_name} ({len(passing)} passing; id_human={ident:.3f})",
            }
        )

    return fasta_lines, provenance_rows, all_audit, updated_df


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional path to config.yaml",
    )
    return parser.parse_args()


def main() -> None:
    """Execute Phase 1 selection, validation, and provenance export."""
    args = parse_args()
    cfg = load_config(args.config)
    log = get_logger("10_select_and_fetch", cfg)
    log.info("Phase 1 goal: fetch and verify somatic cytochrome c for 15 vertebrates + 2 outgroups")

    species_path = REPO_ROOT / cfg["paths"]["species_table"]
    bak_path = species_path.with_suffix(".csv.bak")
    if not bak_path.exists():
        bak_path.write_text(species_path.read_text(encoding="utf-8"), encoding="utf-8")
        log.info("Created backup of species table at %s", bak_path)

    species_df = pd.read_csv(species_path, dtype={"candidate_accession": str})
    raw_dir = resolve(cfg, "raw")
    uniprot_dir = raw_dir / "uniprot"
    uniprot_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = resolve(cfg, "tables")

    human_entry = fetch_uniprot_entry(HUMAN_ANCHOR_ACC, cfg)
    human_seq = human_entry["sequence"]["value"]
    log.info("Fetched human reference anchor %s (length=%d)", HUMAN_ANCHOR_ACC, len(human_seq))

    fasta_lines, prov_rows, audit_rows, updated_df = process_species(
        species_df, cfg, human_seq, uniprot_dir, log
    )

    # Check accession uniqueness across species
    accs = [row["accession"] for row in prov_rows]
    if len(accs) != len(set(accs)):
        raise ValueError(f"Duplicate accessions selected across species: {accs}")

    fasta_path = raw_dir / "cytc_raw.fasta"
    fasta_path.write_text("".join(fasta_lines), encoding="utf-8")

    audit_path = raw_dir / "candidates_all.csv"
    pd.DataFrame(audit_rows).to_csv(audit_path, index=False, encoding="utf-8")

    prov_path = tables_dir / "dataset_provenance.csv"
    pd.DataFrame(prov_rows).to_csv(prov_path, index=False, encoding="utf-8")

    updated_df.to_csv(species_path, index=False, encoding="utf-8")

    log.info(
        "SUMMARY: retrieved %d/%d species (%d candidates evaluated); wrote %s, %s, %s",
        len(prov_rows),
        len(species_df),
        len(audit_rows),
        fasta_path,
        audit_path,
        prov_path,
    )


if __name__ == "__main__":
    main()

