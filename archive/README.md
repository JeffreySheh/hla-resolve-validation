# Archive

The original pipeline and its results, kept as they were. See "History" in the top-level README for what changed and why.

- `scripts/` — `2_merge.sbatch` (chemistry-count merge rule), `3_type.sbatch` (per-sample typing), `grade_1kg.py` (first 1000 Genomes grader), `test_against_truth.py` (first Lai grader; its functions now live in `scripts/test_against_truth.py`), `run_array.sh`/`run_single.sh` (earlier runners)
- `results/` — `hla_validation_final.csv` (Lai, 40 samples), `hla_validation_1kg.csv` (1000 Genomes, 88 samples), `results_summary.csv`
- `data/samples_42.tsv` — the original sample list
- `methods.md` — the original methods write-up
- `matt_download_pilot/` — the download pilot that preceded `1_download.sbatch`
