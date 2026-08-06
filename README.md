# HLA-Resolve Validation

**Independently validating a genomics pipeline that reads human immune-system genes from DNA sequencing benchmarked against a published gold-standard set.**

## What this is

**HLA** genes are the most variable region in the human genome and drive immune recognition — central to transplant matching and autoimmune disease. Reading a person's exact HLA genes from sequencing data ("HLA typing") is notoriously hard.

`hla_resolve` (a pipeline by Matt Glasenapp, Cornejo Lab) does this typing from whole-genome data. **This project independently checks its accuracy** by running it on 39 human samples and comparing every call against a reference set where the true answer is already known (Lai et al., four-field HLA truth).

## Result

**Across 39 samples, the pipeline matched the reference almost perfectly at the protein level. The rare deeper mismatches trace to expected reference-database noise, not typing errors.**

*Concordance = agreement with the reference. Measured at three- and four-field resolution — deeper fields = a stricter test.* Full per-sample results: [`results/hla_validation_results.csv`](results/hla_validation_results.csv).

## Repository layout

```
data/      Lai reference set + sample manifest
scripts/   SLURM jobs to download, merge, and type samples
results/   full results table + readable summary
methods.md detailed methodology and reproduction steps
```

## Credit

Pipeline `hla_resolve` by Matt Glasenapp (Cornejo Lab, UC Santa Cruz). Sequencing data from the Human Pangenome Reference Consortium. Reference HLA types from Lai et al.
