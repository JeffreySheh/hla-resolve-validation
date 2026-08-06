# Methods

How the validation was run, from raw sequencing data to graded HLA calls.

## Data

- **Samples:** 41 of the 44 in the Lai et al. reference set (three were out of scope).
- **Sequencing:** public PacBio HiFi whole-genome data from the Human Pangenome Reference Consortium (HPRC), pulled from S3.
- **Truth set:** Lai et al. four-field HLA types (`data/lai_answer_key.tsv`).

## Workflow (per sample)

1. **Download** the sample's HiFi sequencing runs from the HPRC S3 bucket.
2. **Merge** same-library runs into one BAM to reach enough depth over the HLA region. (Runs from different instruments/libraries are not merged.)
3. **Type** with `hla_resolve` in whole-genome (`WGS`) mode.
4. **Check depth** with `mosdepth` — the mean HLA-region depth is the authoritative coverage number.
5. **Grade** the calls against the Lai reference at three- and four-field resolution.

Batch runs were driven by a sample manifest (`data/samples.tsv`) through a SLURM array job (`scripts/run_array.sh`); single samples used `scripts/run_single.sh`.

## Depth gate

Samples below roughly ⟨10–15⟩× mean HLA-region coverage "soft-fail" — the pipeline types fewer than all 8 genes, and low-depth calls are unreliable. Merging additional runs raises depth and clears the gate.

## Grading

Calls were compared to the Lai reference **by hand**, per gene, at three- and four-field resolution:

- Each gene has two alleles; pairing is unordered (best match wins).
- **Concordance** = how many alleles match the reference.
- Where Lai reports an allele to fewer fields than the call, matching is done at the shallower shared depth (partial credit).

## Interpreting the results

- **Three-field concordance** reflects agreement at the protein-relevant level.
- **Four-field losses** are almost always fourth-field-only: the call is correct through three fields and differs only at the non-coding fourth field — expected reference-database noise, not a typing error.
- A few genuine deeper discrepancies (some DRB1/DQB1 cases) are flagged in the notes column of `results/hla_validation_results.csv`.
