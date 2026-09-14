#!/usr/bin/env python3
"""Reconstruct the PacBio input files behind each 1kGP typing run.

Emulates the two pipeline steps without touching the cluster:

  1_download.sbatch   every non-fail_reads .bam under
                      s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/PacBio_HiFi/
                      for the first cohort (HPRC, then HPRC_PLUS) that has any
  2_merge.sbatch      m84* files if there are at least REVIO_MIN of them,
                      else m64* files, else the m84* files ("Revio(thin)");
                      FORCE_SQ2 samples always take m64*

The result is checked against the read groups hla_resolve wrote to each run
log (data/movies_1kgp.tsv). Where the log disagrees with the rule, the log wins,
because it records what actually ran. Where S3 holds two encodings of one
movie (.5mc / .fixed_kinetics, .hifi_reads / .pm_v1.4.0) the download fetches
both; one filename per movie is listed and the row is flagged. Run logs
confirm only one copy per movie entered the merge.

Library identity comes from HPRC's file index (data/HPRC_PacBio_HiFi.file.index.csv,
github.com/human-pangenomics/HPRC_metadata). Its library_ID carries the movie
name as a suffix, which is stripped so files from one prep share a name.
The index labels every NA19909 file "NA19909_lib1" although the BAM headers
name four preps including two Fiber-seq; the header is right for that sample.

Reads   results/kgp_calls.csv       samples to cover
        data/movies_1kgp.tsv        movie names per sample, from the run logs
        data/HPRC_PacBio_HiFi.file.index.csv
Writes  data/kgp_input_manifest.csv sample, cohort, instrument, library,
                                    file_coverage_x, filename, source,
                                    second_encoding_on_s3
"""

import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CALLS = os.path.join(REPO_ROOT, "results", "kgp_calls.csv")
MOVIES = os.path.join(REPO_ROOT, "data", "movies_1kgp.tsv")
OUT = os.path.join(REPO_ROOT, "data", "kgp_input_manifest.csv")
HPRC_INDEX = os.path.join(REPO_ROOT, "data", "HPRC_PacBio_HiFi.file.index.csv")

S3 = "https://human-pangenomics.s3.amazonaws.com/?list-type=2&prefix=working/%s/%s/raw_data/PacBio_HiFi/"
COHORTS = ("HPRC", "HPRC_PLUS")
REVIO_MIN = 3
FORCE_SQ2 = {"HG01952", "HG01978"}

MOVIE = re.compile(r"(m\d+[a-zA-Z]*_\d+_\d+(?:_s\d)?)")


def s3_listing(sample):
    """(cohort, [filenames]) exactly as 1_download.sbatch would see them."""
    for cohort in COHORTS:
        xml = subprocess.run(["curl", "-s", S3 % (cohort, sample)],
                             capture_output=True, text=True).stdout
        keys = [k.split("/")[-1] for k in re.findall(r"<Key>([^<]+\.bam)</Key>", xml)
                if "fail_reads" not in k]
        if keys:
            return cohort, sorted(keys)
    return "", []


def merge_rule(sample, files):
    revio = [f for f in files if f.startswith("m84")]
    sq2 = [f for f in files if f.startswith("m64")]
    if sample in FORCE_SQ2 and sq2:
        return sq2
    if len(revio) >= REVIO_MIN:
        return revio
    return sq2 or revio


def one_per_movie(files):
    """Prefer the primary encoding when S3 carries a movie twice."""
    best = {}
    for f in files:
        m = MOVIE.match(f).group(1)
        rank = (("fixed_kinetics" in f) + ("pm_v" in f), f)
        if m not in best or rank < best[m][0]:
            best[m] = (rank, f)
    return [best[m][1] for m in sorted(best)], len(files) - len(best)


def library_name(library_id):
    """HPRC appends the movie to library_ID; drop it to get the prep name."""
    return re.sub(r"[._]?m(84|64|54)\d*[a-zA-Z]*_\d+_\d+(_s\d)?.*$", "", library_id)


def main():
    samples = [r["sample"] for r in csv.DictReader(open(CALLS))]
    index = {r["filename"]: r for r in csv.DictReader(open(HPRC_INDEX))}
    log_movies = {}
    for line in open(MOVIES):
        s, _, m = line.rstrip("\n").partition("\t")
        log_movies[s] = set(m.split())

    rows, overrides = [], []
    for s in samples:
        cohort, listing = s3_listing(s)
        chosen = merge_rule(s, listing)
        source = "merge rule"
        logged = log_movies.get(s) or set()
        if logged and logged != {MOVIE.match(f).group(1) for f in chosen}:
            chosen = [f for f in listing if MOVIE.match(f).group(1) in logged]
            source = "run log (differs from merge rule)"
            overrides.append(s)
        elif s in log_movies and not logged and len(chosen) > 1:
            # hla_resolve prints read groups only when the BAM has more than
            # one, so an empty entry means a single input file. The rule's
            # multi-file pick cannot be what ran; the one m84 file is the only
            # single-file candidate, and depth per GB matches it for all nine.
            chosen = [f for f in listing if f.startswith("m84")][:1]
            source = "single read group in run log; sole m84 file"
            overrides.append(s)
        files, n_dup = one_per_movie(chosen)
        prefixes = {f[:3] for f in files}
        inst = ("Revio" if prefixes == {"m84"} else
                "Sequel II" if prefixes == {"m64"} else "mixed")
        for f in files:
            m = index.get(f, {})
            rows.append([s, cohort, inst, library_name(m.get("library_ID", "")),
                         m.get("coverage", ""), f, source, "yes" if n_dup else ""])

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "cohort", "instrument", "library", "file_coverage_x",
                    "filename", "source", "second_encoding_on_s3"])
        w.writerows(rows)

    n = len({r[0] for r in rows})
    print("wrote %s: %d files for %d samples" % (OUT, len(rows), n))
    print("log overrides the merge rule for: %s" % (" ".join(overrides) or "none"))
    libs = defaultdict(set)
    for r in rows:
        if r[3]:
            libs[r[0]].add(r[3])
    multi = sorted(s for s, v in libs.items() if len(v) > 1)
    print("samples merging more than one library prep: %s" % (" ".join(multi) or "none"))
    missing = sorted(set(samples) - {r[0] for r in rows})
    if missing:
        sys.exit("ERROR: no files resolved for %s" % missing)


if __name__ == "__main__":
    main()
