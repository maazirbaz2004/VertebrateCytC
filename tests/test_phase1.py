import re
import pandas as pd
from Bio import SeqIO

from cytc.utils import REPO_ROOT


def test_fasta_matches_provenance():
    recs = list(SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta"))
    prov = pd.read_csv(REPO_ROOT / "results/tables/dataset_provenance.csv")
    assert len(recs) == len(prov) == 17
    assert prov["download_date"].notna().all()
    assert prov["url"].notna().all()
    for acc in prov["accession"]:
        assert (REPO_ROOT / f"data/raw/uniprot/{acc}.json").exists(), acc


def test_required_lineages_present():
    sp = pd.read_csv(REPO_ROOT / "config/species.csv")
    groups = set(sp.loc[sp["candidate_accession"].notna(), "group"])
    for g in ["Primate", "Rodent", "Bird", "Reptile", "Amphibian", "Fish", "Fungus"]:
        assert g in groups, g


def test_every_record_has_heme_motif():
    std_alpha = set("ACDEFGHIKLMNPQRSTVWY")
    for r in SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta"):
        seq = str(r.seq)
        assert re.search(r"C..CH", seq), r.id
        assert set(seq) <= std_alpha, r.id


def test_unique_accessions_and_human_anchor():
    prov = pd.read_csv(REPO_ROOT / "results/tables/dataset_provenance.csv")
    assert prov["accession"].nunique() == len(prov) == 17

    recs = {r.id.split("|")[0]: str(r.seq) for r in SeqIO.parse(REPO_ROOT / "data/raw/cytc_raw.fasta", "fasta")}
    hum = recs["Homo_sapiens"]
    assert len(hum) == 105
    assert hum.startswith("MGDVEKGKKIF")
    assert "CSQCH" in hum
