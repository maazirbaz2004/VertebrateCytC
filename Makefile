PY=python
.PHONY: env check phase1 phase2 phase3 phase4 phase5 phase6 phase7 phase8 all test clean-results

check:
	$(PY) scripts/00_check_env.py

phase1:
	$(PY) scripts/10_select_and_fetch.py
phase2:
	$(PY) scripts/20_prepare_sequences.py
phase3:
	$(PY) scripts/30_align.py
phase4:
	$(PY) scripts/40_conservation.py
phase5:
	$(PY) scripts/50_similarity.py
phase6:
	$(PY) scripts/60_phylogenetics.py
phase7:
	$(PY) scripts/70_interpretation.py
phase8:
	$(PY) scripts/80_make_report_assets.py

all: phase1 phase2 phase3 phase4 phase5 phase6 phase7 phase8

test:
	pytest -q

clean-results:
	rm -rf results/alignments/* results/conservation/* results/similarity/* \
	       results/trees/* results/interpretation/* results/figures/* results/tables/*

