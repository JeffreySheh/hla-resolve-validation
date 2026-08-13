# PacBio HiFi download and merge pilot

A working, verified pipeline for pulling HPRC PacBio HiFi data onto the elkhorn cluster
and combining it into `hla_resolve` inputs. Run end to end on 2026-08-10 for 17 samples
and 1.98 TB.

This document is written to be handed to someone (or some agent) who has to do the same
thing again. It records what was done, what the measured numbers were, what went wrong,
and the cluster-specific facts that are not discoverable from the code.

---

## 1. Result

| stage | outcome |
| --- | --- |
| download | 17 samples, 69 files, 1.98 TB, all byte-exact verified |
| merge | 17 combined bams in `merged/`, all read groups preserved |
| download wall time | ~30 min at ~1.15 GB/s across 17 nodes |
| merge wall time | ~90 min, limited by beegfs at ~800 MB/s combined read and write |

Data lives at `/hb/scratch/mglasena/pacbio_bam/` with one directory per sample holding
the per-movie files, and `merged/<SAMPLE>.hifi.bam` holding the combined inputs.

`hla_resolve` has **not** yet been run on any of it.

---

## 2. The single most important idea

**Downloading is limited by the network card, not the CPU.** Every compute node has the
same 1 GbE external interface, so a 112-core node downloads at exactly the same speed as
a 4-core node. Asking for a bigger machine buys nothing.

The resource you are actually requesting is **machines**. Cores are just the price of
admission. This was measured, not assumed:

- one 112-core node, 12 concurrent files, 192 connections: **100 MB/s**
- 17 nodes, 4 cores each, 4 concurrent files each: **~1150 MB/s**

There are three ceilings stacked on top of each other, and you need all three in mind:

| ceiling | value here | consequence |
| --- | --- | --- |
| per TCP stream | ~45 MB/s | need 3-4 concurrent files to fill one node |
| per node NIC | ~100-115 MB/s (1 GbE) | more cores or more streams do nothing past this |
| site border link | ~1.1-1.2 GB/s (10 GbE) | no point going past ~10-12 saturated nodes |

For CPU-bound work the logic inverts. The merge job asks for 8 cores and does **not** use
`--exclusive`, because there sharing a node is fine and cores are the scarce thing.

---

## 3. Cluster facts (elkhorn, UCSC hummingbird)

These cost hours to discover. They are not in any documentation we found.

**`srun` does not work. Use `salloc`.** Interactive allocations must go through
`salloc -p <part> -A <acct> -q <qos> -w <node> -N1 -c N --mem=NG -t H:MM:SS` and then run
commands in the returned shell. This also rules out `srun --overlap --jobid=...` for
inspecting a node while its job runs, so node diagnosis has to wait for the allocation to
be released.

**Partition sizes.** `windfall` is the only partition with real width, about 88 nodes
(`--partition=windfall --account=windfall --qos=windfall`). It is preemptible, so pair it
with `--requeue`. `lab-colibri` has 2 nodes (node-38, node-39) and `lab-colibri-hmem` has
1 (node-40), both under `--account=pi-jkoc --qos=pi-jkoc`. Nodes are 112 cores.

**Slurm packs array tasks onto one node** if they fit, which silently destroys the entire
point of fanning out for network-bound work. Seventeen 4-core tasks all landed on a single
112-core node and ran no faster than one task. `--exclusive` is what forces one task per
machine. Always verify:

```
squeue -u $USER -o '%.14i %.8T %M %N'
```

and count distinct node names. This is the single easiest mistake to make.

**Nodes 29 and 33 through 36 have 100 Mbit network cards.** Confirmed with
`cat /sys/class/net/*/speed` returning `100` rather than `1000`. They ran at ~11.4 MB/s,
about eight times slower than every other node, and stalled a job that all 12 other nodes
finished. Always exclude them from network-bound work:

```
--exclude=node-29,node-33,node-34,node-35,node-36
```

Known-good 1 GbE nodes from completed transfers: 40, 41, 42, 43, 49, 50, 76, 94, 95, 102,
103, 116, 117. The rest of windfall is unsurveyed. A `-1` in the speed list is InfiniBand
and a `200000` is the 200 Gb IB that carries beegfs, so filesystem traffic is never the
constraint for downloads.

**`aws` lives behind `module load aws`** (2.13.15, at
`/hb/software/apps/aws/prebuilt-2.13.15/bin/aws`). Batch shells do not get the `module`
function for free, so a script has to source the lmod init first.
`download_pacbio.sbatch` does this already.

**samtools comes from the `hla_resolve` conda env**, not a module. `sbatch` propagates
your submitting environment by default, so activating the env before submitting is enough.

**Scratch is `/hb/scratch/mglasena` on beegfs**, 3.7 P with ~2.3 P free. Space is not a
constraint. Note that beegfs is not fast for this workload: the merge measured about
800 MB/s of combined read and write across 17 concurrent tasks.

**`du` under-reports on beegfs** by a few percent because it counts allocated blocks. Use
`du --apparent-size` or a verified file will look incomplete against its expected size.

---

## 4. The files

| file | role |
| --- | --- |
| `make_manifest.sh` | sample list to validated manifest with exact S3 byte sizes |
| `pacbio_manifest.tsv` | the 17-sample manifest actually used (regenerate your own) |
| `fetch_bam.sh` | downloads and verifies one file, called by the driver via `xargs` |
| `download_pacbio.sbatch` | array job, one sample per task |
| `check_downloads.sh` | progress table, per-node rates, orphan cleanup |
| `merge_bams.sbatch` | array job, combines each sample's movies into one bam |

### Paths you must edit

Five lines carry `/hb/scratch/mglasena/pacbio_bam`:

- `download_pacbio.sbatch` line 3 (`#SBATCH --output=`) and line 32 (`OUTDIR`)
- `merge_bams.sbatch` line 3 (`#SBATCH --output=`) and line 43 (`OUTDIR`)
- `check_downloads.sh` line 12 (`OUTDIR`)

The `#SBATCH --output=` lines **cannot** use a variable, so they must be edited directly,
and the logs directory has to exist before you submit because Slurm opens that file
itself. The `OUTDIR` lines are defaults that can alternatively be overridden by
environment variable. You will also need your own `--partition`, `--account`, `--qos`,
and `--array=1-N`.

---

## 5. How to run it

```
module load aws
conda activate hla_resolve          # for samtools
mkdir -p /your/scratch/pacbio_bam/logs
cd /your/scratch/pacbio_bam

./make_manifest.sh samples.tsv > pacbio_manifest.tsv

sbatch --exclusive --exclude=node-29,node-33,node-34,node-35,node-36 download_pacbio.sbatch
./check_downloads.sh -w             # stops itself at N/N

sbatch merge_bams.sbatch            # no --exclusive, this one is CPU and disk bound
STATUS_ONLY=1 bash merge_bams.sbatch
```

`samples.tsv` is tab separated: sample name, cohort (`HPRC` or `HPRC_PLUS`, which is just
the S3 subdirectory), and an optional third column of space-separated movie filenames to
restrict to. Leave the third column off to take every bam in the directory.

```
HG00438	HPRC	m64043_200710_174426.ccs.bam m64043_200711_235708.ccs.bam
HG01109	HPRC_PLUS
```

Both jobs are **safe to rerun**. Verified work is skipped, so a rerun after preemption or
failure only does what is left. This matters on windfall, where preemption is expected.

Useful environment knobs:

| variable | applies to | effect |
| --- | --- | --- |
| `STATUS_ONLY=1` | both | report and exit, changes nothing |
| `SAMPLES="A B"` | both | restrict to some samples |
| `FORCE=1` | download | re-download even if already verified |
| `NJOBS` | download | concurrent files per task, default 4 |
| `VERIFY_MD5=1` | download | check md5 sidecars, costs a full re-read |
| `DELETE_INPUTS=1` | merge | remove per-movie files after each sample verifies |
| `VERIFY_COUNTS=1` | merge | reconcile read counts, roughly doubles runtime |
| `USE_MERGE=1` | merge | force `samtools merge` instead of `cat` |

---

## 6. What gets verified

Nothing is trusted. Each downloaded file must pass all of:

1. **byte-exact size** against the manifest, whose sizes come from S3 itself
2. **`samtools quickcheck -u`**
3. **the 28-byte BGZF end-of-file marker**
   (`1f8b08040000000000ff0600424302001b0003000000000000000000`)

`VERIFY_MD5=1` adds a check against the md5 sidecars S3 publishes next to every object.
It is off by default because aws already checksums each part in flight, and size plus the
EOF marker catches truncation, which is the realistic failure mode.

Each merged file must pass:

1. `quickcheck -u`
2. the BGZF EOF marker
3. **`@RG` count equal to the number of movies** that went into it
4. output size within 3% of the summed inputs

Check 3 is the important one. See section 8.

A file is only considered done when a `.done` marker exists **and** the file it vouches
for is still present. A marker alone is not trusted, so deleting a bam correctly causes a
refetch rather than a false "complete".

---

## 7. Gotchas that cost real time

**`samtools quickcheck` needs `-u` on this data.** These are unaligned bams with no `@SQ`
lines, so plain `quickcheck` returns exit code 8 ("no targets in header") and declares
every single file corrupt. Without `-u` all 69 downloads would have been rejected and
each would have burned three retry attempts before the job gave up. This is the one defect
that would have cost the entire run rather than producing a useful error.

Do not try to work around it by masking bit 8 from the exit code. That was tested and it
wrongly accepts a file that is not a bam at all, because bit 2 gets conflated.

**`aws s3 cp` writes to its own temp name** `<dest>.<random>`, not `<dest>`. Any cleanup
glob has to be `"$f".part*` and not `$f.part`. Getting this wrong orphaned **421 GB** of
abandoned temp files across two cancelled jobs, which nothing would ever have collected.
`check_downloads.sh -c` clears them, and only from samples where every file is verified,
so it can be run safely while a job is still going.

**`aws s3 cp` cannot resume.** A cancelled or preempted transfer loses its partial file
entirely. This is why `NJOBS` matters on a preemptible partition: with `NJOBS=4` a
preemption loses four partial files instead of one. `NJOBS=1` would cap the loss, and
costs nothing in throughput when you are uplink-limited rather than per-stream-limited.

**`du` counts in-flight `.part` files**, so a sample can read near 100% while its files are
still unverified. The verified-file count is the honest progress metric, not bytes.

**Slurm `CG` state is normal.** After `scancel`, tasks sit in COMPLETING while killed
samtools processes flush buffers to beegfs. It clears on its own. If it persists past ten
minutes the nodes are wedged and that is a sysadmin problem.

---

## 8. Why merging uses `samtools cat` and not `samtools merge`

This is the subtlest thing in the pipeline and the easiest to get silently wrong.

Each PacBio movie file has its own `@RG` (read group) line with a distinct ID, and every
read carries an `RG:Z:` tag pointing at it. When you combine movies you must preserve all
of them.

**`samtools merge` is correct but slow.** It handles headers properly, but it decodes and
re-encodes every record. Measured here at **15 MB/s per task**, which projected to a 4 to
5 hour job for 2 TB.

**Bare `samtools cat` is fast but silently wrong.** It copies compressed BGZF blocks
untouched, so it runs at I/O speed. But it takes the header from the **first file only**,
which drops the `@RG` lines of every later movie and leaves those reads pointing at read
groups that do not exist in the header. Demonstrated on real HG00438 data:

```
samtools cat   -> 1 @RG in header, but 2 distinct RG:Z tags on reads   # broken
samtools merge -> 2 @RG in header, 2 distinct RG:Z tags on reads       # correct
```

**The fix is `samtools cat -h <header>`** with a header you build yourself containing every
movie's `@RG`. `merge_bams.sbatch` constructs it as `@HD` from the first file, then every
distinct `@RG` across all movies, then the remaining header lines from the first file.

This was verified equivalent to `samtools merge` on real data: identical read names,
sequences, qualities, and read-to-read-group mapping. The only difference is tag ordering
within each record, which carries no meaning in BAM.

The script refuses to use the fast path if the constructed header does not end up with
exactly one `@RG` per movie, falling back to `samtools merge` automatically. It also falls
back if kinetics stripping is needed, since that requires decoding every record anyway.

**Expect ~99% of summed input bytes** in the output. Concatenation drops the duplicate
headers and EOF blocks of all but one movie, so the result is legitimately a hair smaller.

---

## 9. Kinetics

PacBio kinetics tags (`fi`, `fp`, `ri`, `rp`, `ipd`, `pw`) roughly double file size and
nothing downstream uses them.

**The HPRC raw `ccs.bam` release has them already stripped.** Six files were sampled across
both cohorts and every read carries only `RG np rq sn zm`, with one batch adding
`ec we ws`. So there is nothing to strip and no space to save. What actually makes these
files 30 GB apiece is per-base quality strings on ~19 kb reads, and those must be kept
because DeepVariant uses them.

`merge_bams.sbatch` checks anyway and will strip them inside the merge pipe if it ever
finds them, at the cost of forcing the slow decode path.

---

## 10. Data source, and a warning about other samples

Files come from `s3://human-pangenomics/working/{HPRC|HPRC_PLUS}/{sample}/raw_data/PacBio_HiFi/`.
Everything is also reachable over plain HTTPS at
`https://human-pangenomics.s3.amazonaws.com/...`, which is useful for quick checks with
`curl` when `aws` is not loaded. md5 sidecars sit next to every object.

**The 17 samples here are uniformly Sequel II `m64xxx_*.ccs.bam`. Newer samples are not.**
Checked on 2026-08-10, other HPRC sample directories mix three things:

- `m64xxx_*.ccs.bam` — Sequel II, what this pilot used
- `m84xxx_*.hifi_reads.bcNNNN.bam` — Revio, newer chemistry, often much larger files
- `m84xxx_*.fail_reads.bcNNNN.bam` — **QC-failed reads, must be excluded**

`make_manifest.sh` in its current form takes every `.bam` in a directory, so **it needs a
`fail_reads` filter before being pointed at those samples**. It would also silently mix two
chemistries in one merged file, which may or may not be what you want.

For a concrete example, an 18-sample list checked on 2026-08-10 (HG02109, HG02055,
HG01928, HG02145, HG02723, HG02257, HG02572, HG01891, HG03516, HG01361, HG01071, HG01952,
HG01978, HG01175, HG01358, HG02148, HG00741, HG00673) broke down as 1.82 TB of Sequel II,
3.37 TB of Revio, and 0.25 TB of `fail_reads`. HG01361 has Revio only and no Sequel II at
all, so a Sequel-only cohort from that list would be 17 samples rather than 18.

---

## 11. Diagnosing a slow or stuck job

Ordered by how often each turned out to be the answer.

**Are the tasks actually on different nodes?**

```
squeue -u $USER -o '%.14i %.8T %M %N'
```

Count distinct node names. If they are stacked, you forgot `--exclusive`.

**Is one node slower than the others?**

```
./check_downloads.sh -r
```

Prints MB/s per node over a 60-second window. ~100 MB/s is a saturated 1 GbE and is what
you want. ~11 MB/s means a 100 Mbit node: cancel that task and resubmit with the node
added to `--exclude`.

**Is it progressing at all?** For the merge, watch the `.tmp` files grow:

```
cd merged
stat -c '%n %s' *.tmp | sort > /tmp/m1; sleep 30
stat -c '%n %s' *.tmp | sort > /tmp/m2
join /tmp/m1 /tmp/m2 | awk '{printf "%-26s %6.1f MB/s\n", $1, ($3-$2)/30/1000000}'
```

**What phase is each task in?** The last log line tells you:

```
for i in $(seq 1 17); do
  printf '%2d  %s\n' $i "$(tail -n 1 logs/merge_<JOBID>_${i}.out)"
done
```

**Have any files failed verification?**

```
grep -h 'FAILED\|SKIP\|WARNING\|attempt [23]' logs/*.out
```

Silence means everything is on its first attempt and verifying clean.

A cautionary note on method. During this pilot a conclusion was drawn that five stalled
tasks had written zero bytes, based on comparing a stale `du` reading against manifest
totals that happened to match by coincidence. The actual directory listing showed they
were downloading, just slowly. **Look at the files before theorizing about them.**

---

## 12. Final verification

After the merge finishes, from the pacbio_bam directory with samtools available:

```
echo "outputs: $(ls merged/*.hifi.bam 2>/dev/null | wc -l)"
echo "leftover tmp: $(ls merged/*.tmp 2>/dev/null | wc -l)"
grep -h 'FAILED\|SKIP\|WARNING' logs/merge_*.out
samtools quickcheck -u merged/*.hifi.bam && echo "all structurally valid"

for f in merged/*.hifi.bam; do
  s=$(basename "$f" .hifi.bam)
  want=$(awk -F'\t' -v x="$s" '$1==x' pacbio_manifest.tsv | wc -l)
  tgt=$(awk -F'\t' -v x="$s" '$1==x{t+=$4} END{print t}' pacbio_manifest.tsv)
  got=$(samtools view -H "$f" | grep -c '^@RG')
  cur=$(stat -c %s "$f")
  printf '%-9s %s movies %s @RG  %3d%%  %s\n' "$s" "$want" "$got" $((cur*100/tgt)) \
    "$([ "$want" = "$got" ] && echo OK || echo RG-MISMATCH)"
done
```

You want every sample at `N movies N @RG 99% OK`. The `@RG` column is the one that would
have caught a bare-`cat` mistake, so it is worth reading rather than skimming.

Do not run `DELETE_INPUTS=1` until `hla_resolve` has successfully consumed at least one
merged file. Re-downloading is half an hour, and keeping the inputs costs nothing against
2.3 PB of free scratch.

---

## 13. Sample set used

17 samples, 69 movie files, 1.98 TB. Eight HPRC and nine HPRC_PLUS, which is only an S3
path difference and not a data difference.

```
HG00438  HG00735  HG01106  HG01109  HG01243  HG02622  HG02630  HG02717  HG02818
HG02886  HG03098  HG03486  HG03540  NA18906  NA19240  NA20129  NA21309
```

Most samples have 4 movie files. HG02818 has 7 and NA20129 has 5. HG01109, HG01243, and
HG03098 have 3. Per-sample volume ranges from 94.8 GB (NA20129) to 139.2 GB (HG01106).
