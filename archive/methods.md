# Methods

How the validation was run, from raw sequencing data to graded HLA calls.

## Data

- **Samples:** 41 of the 44 in the Lai et al.
- **Sequencing:** public PacBio HiFi whole-genome data from the Human Pangenome Reference Consortium (HPRC), pulled from S3.
- **Truth set:** Lai et al. four-field HLA types (`data/lai_answer_key.tsv`).

## Workflow (per sample)

1. **Download** the sample's HiFi sequencing runs from the HPRC S3 bucket.
2. **Merge** same-library runs into one BAM to reach enough depth over the HLA region. (Runs from different instruments/libraries are not merged.)
3. **Type** with `hla_resolve` in whole-genome (`WGS`) mode.
4. **Check depth** with `mosdepth` the mean HLA-region depth is the authoritative coverage number.

Batch runs were driven by a sample manifest (`data/samples.tsv`) through a SLURM array job (`scripts/run_array.sh`); single samples used `scripts/run_single.sh`.

## Depth gate

hla_resolve applies a per-gene depth gate: a gene is typed only when its antigen recognition site (ARS) reaches ≥ 8× mean coverage depth. Genes below that threshold soft-fail, and their calls are unreliable. Merging additional same-library runs raises per-gene ARS depth and clears the gate.

## Grading

Calls were compared to the Lai reference **by hand**, per gene, at three- and four-field resolution:

- Each gene has two alleles; pairing is unordered (best match wins).
- **Concordance** = how many alleles match the reference.
- Where Lai reports an allele to fewer fields than the call, matching is done at the shallower shared depth (partial credit).

## Interpreting the results

- **Three-field concordance** requires agreement not only at the protein level (two fields — the amino acid sequence) but also on synonymous coding substitutions (the third field — silent DNA changes that don't alter the protein). It is therefore a stricter bar than protein-level identity.
- **Four-field losses** are almost always fourth-field-only: the call agrees through three fields and differs only at the non-coding fourth field. The three-field call is therefore correct; the fourth-field difference is consistent with reference-database noise at non-coding positions, though we don't rule out a genuine fourth-field difference.
