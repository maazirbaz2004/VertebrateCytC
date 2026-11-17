# PHASE 1: Species Selection and Sequence Retrieval

**Read `CLAUDE.md` first. Prerequisite: Phase 0 is done.**

## 0. Objective

Assemble the raw dataset: **15 vertebrate species + 2 non-vertebrate outgroups** (17 total), each with one verified somatic cytochrome c protein sequence downloaded from UniProt (primary) or NCBI (fallback), with full provenance recorded.

**Serves:** Objective 1; foundation for RQ1-RQ4.

## 1. Prerequisites

- Phase 0 validation checks pass.
- Internet access to `rest.uniprot.org` and (fallback) `eutils.ncbi.nlm.nih.gov`.
- `config.yaml` has a real `ncbi_email`.

## 2. Inputs

- `config/config.yaml` (`data_collection` section).
- The candidate species table in section 3.1 (to be written to `config/species.csv`).

## 3. Tasks

### Task 1.1: Write the species table `config/species.csv`

Design rationale: the proposal asks for 10-15 vertebrates "from different evolutionary groups including mammals, birds, and other vertebrates". We use **15 vertebrates spanning all major lineages** and **2 outgroups** so the tree can be rooted and divergence is measured over a deep timescale. Eight mammals are included on purpose, so that RQ3 can compare *closely related* pairs (human-macaque, mouse-rat) against *distant* pairs.

Columns: `species, common_name, taxon_id, group, role, approx_mya_from_human, uniprot_gene_hint, candidate_accession, notes`

```csv
species,common_name,taxon_id,group,role,approx_mya_from_human,uniprot_gene_hint,candidate_accession,notes
Homo_sapiens,Human,9606,Primate,ingroup,0,CYCS,P99999,reference species for numbering
Macaca_mulatta,Rhesus macaque,9544,Primate,ingroup,29,CYCS,,closest relative to human in set
Mus_musculus,House mouse,10090,Rodent,ingroup,90,Cycs,P62897,also has testis-specific Cyct: EXCLUDE paralog
Rattus_norvegicus,Brown rat,10116,Rodent,ingroup,90,Cycs,P62898,also has testis-specific paralog: EXCLUDE
Canis_lupus_familiaris,Dog,9615,Carnivore,ingroup,94,CYCS,,
Bos_taurus,Cattle,9913,Artiodactyl,ingroup,94,CYCS,P62894,classic cytochrome c species
Equus_caballus,Horse,9796,Perissodactyl,ingroup,94,CYCS,P00004,historic reference for classical numbering
Monodelphis_domestica,Gray short-tailed opossum,13616,Marsupial,ingroup,160,CYCS,,most divergent mammal here
Gallus_gallus,Chicken,9031,Bird,ingroup,320,CYCS,P67881,
Taeniopygia_guttata,Zebra finch,59729,Bird,ingroup,320,CYCS,,second bird (tests bird monophyly)
Anolis_carolinensis,Green anole,28377,Reptile,ingroup,320,CYCS,,squamate reptile
Chrysemys_picta,Painted turtle,8478,Reptile,ingroup,320,CYCS,,turtle (sister to birds+crocs molecularly)
Xenopus_tropicalis,Western clawed frog,8364,Amphibian,ingroup,352,cycs,,
Danio_rerio,Zebrafish,7955,Fish,ingroup,430,cycs,,teleost; check for multiple cycs paralogs
Takifugu_rubripes,Torafugu,31033,Fish,ingroup,430,cycs,,teleost, second fish
Drosophila_melanogaster,Fruit fly,7227,Insect,outgroup,780,Cyt-c-p,,several paralogs: choose the main somatic one (Cyt-c-p)
Saccharomyces_cerevisiae,Baker's yeast,559292,Fungus,outgroup,1100,CYC1,P00044,iso-1 cytochrome c; N-terminal extension
```

**IMPORTANT.**
- `candidate_accession` values are **from memory and must be verified** by the retrieval script (organism matches, length in range, contains the heme motif, protein name contains "cytochrome c"). Where blank, the script discovers it.
- `approx_mya_from_human` values are rough literature-style divergence times used only for the Phase 5/7 sequence-identity-vs-time analysis. In Phase 7, **replace them with TimeTree (timetree.org) estimates if accessible, and note the source**. Otherwise keep and label as "approximate".
- If a species cannot yield a trustworthy sequence, substitute a close relative (e.g. *Alligator mississippiensis* 8496 for the turtle, *Pan troglodytes* 9598 for macaque) and log the change in `PROGRESS.md` -> Deviations. Keep >= 10 vertebrates, ideally 15, and keep at least one fish, one amphibian, one reptile and one bird.

### Task 1.2: Write `scripts/10_select_and_fetch.py`

Strategy per species (in this order, stop at first success):

1. **Direct accession lookup** if `candidate_accession` is set: fetch `https://rest.uniprot.org/uniprotkb/<ACC>.json` and validate.
2. **UniProt search** by gene + taxon:
   `query=(gene_exact:<GENE>) AND (organism_id:<TAXID>)`
   with `fields=accession,id,gene_names,protein_name,organism_name,organism_id,length,reviewed,sequence`.
   Try `reviewed:true` first, then unreviewed (TrEMBL).
3. **UniProt search by protein name**: `(protein_name:"cytochrome c") AND (organism_id:<TAXID>)`, filtered by length range.
4. **NCBI Protein fallback** via Biopython Entrez: `esearch(db="protein", term='cytochrome c[Title] AND txid<TAXID>[Organism:exp]')` then `efetch(rettype="fasta")`.

**Validation filter applied to every candidate** (all must pass):

| Check | Rule |
|-------|------|
| Organism | `organism_id` equals the species `taxon_id` (for yeast, accept 559292 or 4932 / strain variants but record which) |
| Length | within `[vertebrate_len_min, vertebrate_len_max]` for ingroup, `[outgroup_len_min, outgroup_len_max]` for outgroup |
| Heme motif | regex `C..CH` (from config) found in sequence |
| Name | protein name matches `cytochrome c` (case-insensitive) and does **not** contain `oxidase`, `c1`, `peroxidase`, `reductase`, `-type`, `like` or `testis` |
| Alphabet | only the 20 standard amino acids (no `X`, `B`, `Z`, `U`, `J`) |
| Not a fragment | UniProt `fragment` flag absent (or length OK) |

**Ranking when multiple candidates pass:** (1) reviewed (Swiss-Prot) over unreviewed, (2) exact gene-name match, (3) length closest to 104/105 for vertebrates, (4) lowest accession alphabetically as a deterministic tiebreak. **Log all passing candidates** and which one was picked and why.

**Paralog handling (important).**
- Mouse and rat have a testis-specific cytochrome c (*Cyct*). Exclude it (name contains "testis").
- *Drosophila* has *Cyt-c-p* (the main one, CG17903), *Cyt-c-d* (CG13263), and others. Choose the **Cyt-c-p** gene product; if the search returns several, prefer the one whose gene name matches the hint.
- Zebrafish/Takifugu may have more than one *cycs*-like gene. If two pass, pick the one with the higher identity to the human sequence **only as a last-resort tiebreak** and log it. Better: prefer the gene named exactly `cycs`/`cycsa`.
- Never include pseudogene-derived or "cytochrome c, somatic-like/predicted LOW QUALITY" records if a better one exists.

**Outputs written by the script:**

| File | Content |
|------|---------|
| `data/raw/uniprot/<ACC>.json` | raw API response for each chosen record (untouched) |
| `data/raw/cytc_raw.fasta` | all chosen sequences, header `>Genus_species|ACCESSION|source` |
| `data/raw/candidates_all.csv` | **every** candidate examined (species, accession, passed?, failure reasons) |
| `config/species.csv` | **updated** with the verified `candidate_accession` filled in (keep a backup `species.csv.bak`) |
| `results/tables/dataset_provenance.csv` | species, taxon_id, accession, database (UniProt/NCBI), reviewed flag, gene name, protein name, length, download date (ISO), URL used, notes |

Reference skeleton:

```python
"""Phase 1: select and fetch one cytochrome c sequence per species."""
import json, re, sys, time
from datetime import date
from pathlib import Path

import pandas as pd
import requests
from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cytc.utils import REPO_ROOT, get_logger, load_config, resolve

FIELDS = "accession,id,gene_names,protein_name,organism_name,organism_id,length,reviewed,sequence,fragment"
BAD_NAME = re.compile(r"oxidase|c1\b|peroxidase|reductase|-type|like|testis|subunit", re.I)
STD = set("ACDEFGHIKLMNPQRSTVWY")


def uniprot_search(query: str, cfg: dict, size: int = 25) -> list[dict]:
    r = requests.get(cfg["data_collection"]["uniprot_rest"],
                     params={"query": query, "format": "json", "fields": FIELDS, "size": size},
                     timeout=60)
    r.raise_for_status()
    time.sleep(cfg["data_collection"]["request_pause_sec"])
    return r.json().get("results", [])


def validate(entry: dict, sp: pd.Series, cfg: dict) -> list[str]:
    """Return a list of failure reasons (empty list == passes)."""
    dc, why = cfg["data_collection"], []
    seq = entry["sequence"]["value"]
    pname = (entry.get("proteinDescription", {}).get("recommendedName", {})
             .get("fullName", {}).get("value", ""))
    lo, hi = ((dc["outgroup_len_min"], dc["outgroup_len_max"]) if sp["role"] == "outgroup"
              else (dc["vertebrate_len_min"], dc["vertebrate_len_max"]))
    # NOTE: the JSON field layout of UniProt entries is nested. Inspect one real response
    # first and adapt these accessors (organism taxonId, fragment flag, etc.).
    if not (lo <= len(seq) <= hi): why.append(f"length {len(seq)} outside [{lo},{hi}]")
    if not re.search(dc["heme_motif_regex"], seq): why.append("no CxxCH motif")
    if set(seq) - STD: why.append("non-standard residues")
    if pname and BAD_NAME.search(pname): why.append(f"bad protein name: {pname}")
    if pname and "cytochrome c" not in pname.lower(): why.append(f"name lacks 'cytochrome c': {pname}")
    tax = entry.get("organism", {}).get("taxonId")
    if tax is not None and int(tax) != int(sp["taxon_id"]) and not (sp["species"] == "Saccharomyces_cerevisiae" and int(tax) in (559292, 4932, 4932)):
        why.append(f"taxon mismatch {tax} != {sp['taxon_id']}")
    return why

# ... implement: fetch_by_accession(), search strategies 2-3, ncbi_fallback(),
#     rank_candidates(), main() writing every output file listed above.
```

Implement the full script, **inspect one real UniProt JSON response first** (`curl "https://rest.uniprot.org/uniprotkb/P99999.json" | python -m json.tool | head -150`) to confirm the field paths (`proteinDescription.recommendedName.fullName.value`, `organism.taxonId`, `sequence.value`, `genes[].geneName.value`, `entryType`) before writing accessors.

### Task 1.3: Manual sanity review (the agent does this explicitly)

After the script runs, **print and read** `candidates_all.csv` and `dataset_provenance.csv`. For each species confirm:
- the organism name in the record matches the species;
- the protein is the **main/somatic** cytochrome c, not a paralog;
- the length is plausible (vertebrates about 104-105, yeast CYC1 about 108, *Drosophila* about 105-108).

Known-good anchors for **sanity checking only** (do not hard-code them into the pipeline logic): human CYCS (P99999) is 105 aa, begins `MGDVEKGKKIF` and contains the motif `CSQCH`. If your downloaded human sequence disagrees, something is wrong with the retrieval code.

### Task 1.4: Write the tests `tests/test_phase1.py`

```python
import pandas as pd
from Bio import SeqIO
from cytc.utils import REPO_ROOT

def test_fasta_matches_provenance():
    recs = list(SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta"))
    prov = pd.read_csv(REPO_ROOT / "results/tables/dataset_provenance.csv")
    assert len(recs) == len(prov) >= 12

def test_required_lineages_present():
    sp = pd.read_csv(REPO_ROOT / "config/species.csv")
    groups = set(sp.loc[sp.candidate_accession.notna(), "group"])
    for g in ["Primate", "Rodent", "Bird", "Reptile", "Amphibian", "Fish", "Fungus"]:
        assert g in groups

def test_every_record_has_heme_motif():
    import re
    for r in SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta"):
        assert re.search("C..CH", str(r.seq)), r.id
```

## 4. Expected outputs

| Path | Description |
|------|-------------|
| `data/raw/cytc_raw.fasta` | 15-17 verified sequences |
| `data/raw/uniprot/*.json` | raw API responses |
| `data/raw/candidates_all.csv` | audit trail of all candidates considered |
| `config/species.csv` | accessions filled in |
| `results/tables/dataset_provenance.csv` | provenance table (goes into the report's methods) |
| `logs/10_select_and_fetch_*.log` | run log |

## 5. Validation checks

- [ ] >= 10 vertebrates (target 15) and both outgroups present. Any shortfall is explained in `PROGRESS.md`.
- [ ] At least one each of: primate, rodent, bird, reptile, amphibian, fish, fungus.
- [ ] Every sequence passes the filter table in Task 1.2.
- [ ] No two species share an accession (except by documented exception).
- [ ] Human sequence is 105 aa and matches the sanity anchors.
- [ ] Provenance CSV has a download date and URL for every record.
- [ ] `pytest -q tests/test_phase1.py` passes.

## 6. Common pitfalls

- **UniProt query syntax.** Use `gene_exact:` and `organism_id:` fields, and URL-encode spaces and quotes (let `requests` do it via `params=`).
- **Paralog contamination** (testis-specific, *Drosophila* paralogs): the single most likely source of silent error. A wrong paralog will distort the tree. Check the protein names.
- **Yeast taxon ID.** UniProt files *S. cerevisiae* S288C under taxon 559292. Searching 4932 may return different entries.
- **Predicted/low-quality records** (e.g. "LOW QUALITY PROTEIN") can have truncations. The length and motif filters catch most.
- **Do not trust memory for sequences.** Only accession *candidates* are given above, never residues.
- **Rate limits.** Keep the pause between requests, and use a descriptive `User-Agent` header if you get 429s.
- **Network blocked:** tell the user which domain failed rather than inventing data.

## 7. Done when

- [ ] `make phase1` reproduces all outputs from scratch.
- [ ] Section 5 fully ticked; deviations logged in `PROGRESS.md`.
- [ ] Phase 1 committed.
