#!/usr/bin/env python3
"""Choose the input BAMs to merge for each sample: one library prep, WGS only.

Replaces the chemistry rule in 2_merge.sbatch (m84* count, else m64*) with a
decision made from metadata and written to a plan the merge step consumes.

Inputs
  data/library_by_file.tsv              every downloaded BAM on the cluster:
                                        sample, filename, bytes, and the LB, PU
                                        and SM tags from the BAM's own header
  data/HPRC_PacBio_HiFi.file.index.csv  HPRC file index: library_ID, coverage
                                        (github.com/human-pangenomics/HPRC_metadata)
  data/merge_overrides.tsv              hand-written corrections, see below

Rules, per sample
  1. drop files under MIN_FILE_GB -- barcode cross-talk stubs
  3. drop files whose LB names another sample -- cross-talk that is not tiny
  4. drop non-WGS preps: LB matching fiber|ampure|dilut|test
  5. keep one encoding per movie (.5mc / .fixed_kinetics / .pm_v1.4.0 are the
     same reads); the plain one wins
  6. group the rest by prep: HPRC library_ID with the movie suffix stripped,
     else the header LB, else "<sample>:untagged" (older ccs.bam carry no tags)
  7. rank preps by depth: summed HPRC-index coverage when every file of the
     prep is indexed, else bytes scaled per file type (ccs.bam without
     kinetics tags is ~3x denser than hifi_reads/5mc with them). The deepest
     prep is the sample's primary entry; if it is below MIN_DEPTH the next
     prep is merged in and the entry flagged. Every other prep that reaches
     MIN_DEPTH on its own becomes a further entry: a different library of the
     same sample is an independent measurement and is benchmarked separately.
     Entry ids are <sample>__<prep>, with the prep name made filesystem-safe.
     A prep that spans two instruments is flagged.
  8. cap depth: cells of the chosen prep are added in chronological order
     (movie name) until estimated depth reaches TARGET_DEPTH; once above
     MIN_DEPTH, a cell that would push past MAX_DEPTH is left out; if the
     total still exceeds MAX_DEPTH, the smallest cells are shed while the
     sample stays above the floor. Deep samples otherwise inflate the
     benchmark. Depth here is the HPRC-index
     genome-wide value per file (or the bytes estimate), not measured.

Scope
  select_inputs.py [samples.txt] [--refresh]
  With a sample list (one ID per line) only those samples are planned;
  otherwise every sample in library_by_file.tsv.

Overrides (data/merge_overrides.tsv, tab-separated: sample, action, pattern, note)
  exclude  <regex on filename>            drop matching files before rule 5
  prep     <regex on filename>=<name>     assign a prep name (metadata is wrong)
  depth    <prep name>=<x>                 pin a prep's total depth when neither the
                                           index nor the byte estimate is right

Outputs
  data/merge_plan.tsv           input_bam, sample, prep, filename, size_gb, cov_x, note
  data/merge_plan_excluded.tsv  sample, filename, reason
  --refresh also lists S3 and writes data/download_manifest.tsv in the layout
  1_download.sbatch reads (sample, filename, s3 uri, bytes) for files that are
  on S3 but not on disk, after the filename-level rules. The listing is
  recursive: HPRC keeps some cells in subfolders (primrose/) that a plain
  "aws s3 ls" never shows, and the URI carries the full key. Header-level rules
  (3, 4) can only run once a file is downloaded.
"""

import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = lambda *p: os.path.join(REPO_ROOT, "data", *p)  # noqa: E731

FILES = D("library_by_file.tsv")
INDEX = D("HPRC_PacBio_HiFi.file.index.csv")
OVERRIDES = D("merge_overrides.tsv")
PLAN = D("merge_plan.tsv")
EXCLUDED = D("merge_plan_excluded.tsv")
DOWNLOAD = D("download_manifest.tsv")

MIN_FILE_GB = 1.0       # smaller than any real cell; catches demultiplex spill
MIN_DEPTH = 20.0        # add a second prep only below this
TARGET_DEPTH = 30.0     # stop adding cells once estimated depth reaches this
MAX_DEPTH = 40.0        # never add a cell that would exceed this (above MIN_DEPTH)
# x per GB when a file is not in the HPRC index, measured on this cohort:
# ccs.bam (no kinetics tags) ~0.23; hifi_reads/5mc with kinetics ~0.08
X_PER_GB = {"ccs": 0.23, "default": 0.08}


def est_cov(e):
    if e["cov"] is not None:
        return e["cov"]
    rate = X_PER_GB["ccs"] if e["name"].endswith(".ccs.bam") else X_PER_GB["default"]
    return e["size"] / 1e9 * rate


def instrument(name):
    return "Revio" if name.startswith("m84") else "Sequel II"


def entry_id(sample, preps):
    tag = "+".join(re.sub(r"[^A-Za-z0-9.]+", "_", p).strip("_") for p in preps)
    return "%s__%s" % (sample, tag)
S3 = "https://human-pangenomics.s3.amazonaws.com/?list-type=2&prefix=working/%s/%s/raw_data/PacBio_HiFi/"
S3_URI = "s3://human-pangenomics/%s"      # full key: files may sit in subfolders (primrose/)
COHORTS = ("HPRC", "HPRC_PLUS")

MOVIE = re.compile(r"(m\d+[a-zA-Z]*_\d+_\d+(?:_s\d)?)")
NON_WGS = re.compile(r"fiber|ampure|dilut|test", re.I)
SAMPLE_ID = re.compile(r"(?<![0-9])(?:HG|NA|PG)?(\d{4,5})(?![0-9])")


def movie(f):
    m = MOVIE.match(f)
    return m.group(1) if m else f


def encoding_rank(f):
    """Lower is preferred when one movie has several encodings."""
    return (("fixed_kinetics" in f) + ("pm_v" in f), len(f))


def strip_movie(library_id):
    return re.sub(r"[._]?m(84|64|54)\d*[a-zA-Z]*_\d+_\d+(_s\d)?.*$", "", library_id)


def foreign_library(sample, lb, sm):
    """True if LB or SM names other samples and never this one.
    Pool names list several samples, with or without the letter prefix and
    with or without the leading zero ("HG01978_1928_2148_...")."""
    own = {sample[2:], sample[2:].lstrip("0")}
    for tag in (lb, sm):
        ids = {i.lstrip("0") for i in SAMPLE_ID.findall(tag)}
        if ids and not (own & ids):
            return True
    return False


def load_files():
    """library_by_file.tsv: sample, filename, bytes, LB, PU, SM (header line)."""
    files = defaultdict(list)
    for line in open(FILES):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or not f[2].isdigit():
            continue
        files[f[0]].append(dict(name=f[1], size=int(f[2]), lb=f[3], pu=f[4], sm=f[5]))
    return files


def load_overrides():
    ov = defaultdict(list)
    if not os.path.exists(OVERRIDES):
        return ov
    for line in open(OVERRIDES):
        if not line.strip() or line.startswith("#"):
            continue
        sample, action, pattern, *note = line.rstrip("\n").split("\t")
        ov[sample].append((action, pattern, note[0] if note else ""))
    return ov


def plan_sample(sample, entries, index, overrides):
    kept, dropped = [], []

    def drop(e, why):
        dropped.append((sample, e["name"], why))

    for e in entries:
        if e["size"] < MIN_FILE_GB * 1e9:
            drop(e, "under %g GB: cross-talk stub" % MIN_FILE_GB); continue
        if foreign_library(sample, e["lb"], e["sm"]):
            drop(e, "LB names another sample: %s" % e["lb"]); continue
        if NON_WGS.search(e["lb"]):
            drop(e, "non-WGS prep: %s" % e["lb"]); continue
        hit = next((n for a, p, n in overrides.get(sample, [])
                    if a == "exclude" and re.search(p, e["name"])), None)
        if hit is not None:
            drop(e, "override: %s" % hit); continue
        kept.append(e)

    # one encoding per movie
    by_movie = defaultdict(list)
    for e in kept:
        by_movie[movie(e["name"])].append(e)
    kept = []
    for m, es in by_movie.items():
        es.sort(key=lambda e: encoding_rank(e["name"]))
        kept.append(es[0])
        for e in es[1:]:
            drop(e, "same movie as %s" % es[0]["name"])

    # prep per file
    for e in kept:
        e["cov"] = None
        # the index lists one encoding per movie (usually the .ccs.bam or
        # .hifi_reads.bam); a re-encoded copy of the same movie shares its
        # library and coverage
        meta = index.get(e["name"]) or index.get((sample, movie(e["name"])))
        if meta and meta.get("library_ID"):
            e["prep"] = strip_movie(meta["library_ID"])
            e["cov"] = float(meta["coverage"]) if meta.get("coverage") else None
            e["src"] = "index"
        elif e["lb"]:
            e["prep"] = e["lb"]
            e["src"] = "LB"
        else:
            e["prep"] = sample + ":untagged"
            e["src"] = "untagged"
        for a, p, n in overrides.get(sample, []):
            if a == "prep" and "=" in p:
                rx, name = p.split("=", 1)
                if re.search(rx, e["name"]):
                    e["prep"] = name
                    e["src"] = "override"

    preps = defaultdict(list)
    for e in kept:
        preps[e["prep"]].append(e)
    # pinned depths: spread evenly over the prep's files so the cap still works
    for a, p, n in overrides.get(sample, []):
        if a == "depth" and "=" in p:
            name, x = p.split("=", 1)
            if name in preps:
                for e in preps[name]:
                    e["cov"] = float(x) / len(preps[name])
                    e["src"] = "override (depth pinned)"
    if not preps:
        return [], dropped, "no usable files"

    def depth(p):
        return sum(est_cov(e) for e in preps[p])

    def indexed(p):
        return all(e["cov"] is not None for e in preps[p])

    order = sorted(preps, key=depth, reverse=True)

    # primary entry: deepest prep, plus the next one only if below the floor
    entries = []                      # list of (prep list, notes)
    primary, pnotes = [order[0]], []
    if not indexed(order[0]):
        pnotes.append("depth estimated from bytes (%.0fx)" % depth(order[0]))
    if depth(order[0]) < MIN_DEPTH:
        if len(order) > 1:
            primary.append(order[1])
            pnotes.append("below %gx (%.1fx); second prep %s merged in"
                          % (MIN_DEPTH, depth(order[0]), order[1]))
        else:
            pnotes.append("below %gx (%.1fx); no other prep on disk" % (MIN_DEPTH, depth(order[0])))
    entries.append((primary, pnotes))
    # further entries: any remaining prep that stands on its own
    for p in order:
        if p in primary:
            continue
        if depth(p) >= MIN_DEPTH:
            entries.append(([p], ["additional library"] +
                            ([] if indexed(p) else ["depth estimated from bytes (%.0fx)" % depth(p)])))
        else:
            for e in preps[p]:
                drop(e, "prep %s is %.0fx, under %g on its own" % (p, depth(p), MIN_DEPTH))

    rows, notes_out = [], []
    for chosen, notes in entries:
        cells = [e for p in chosen for e in preps[p]]
        cells.sort(key=lambda e: (chosen.index(e["prep"]), movie(e["name"])))
        insts = {instrument(e["name"]) for e in cells}
        if len(insts) > 1:
            notes.append("prep spans " + " + ".join(sorted(insts)))
        # rule 8: depth cap, chronological within the entry
        taken, total = [], 0.0
        for e in cells:
            c = est_cov(e)
            if taken and total >= TARGET_DEPTH:
                drop(e, "depth cap: already %.0fx of %g target" % (total, TARGET_DEPTH)); continue
            if taken and total >= MIN_DEPTH and total + c > MAX_DEPTH:
                drop(e, "depth cap: %.0fx + %.0fx would exceed %g" % (total, c, MAX_DEPTH)); continue
            taken.append(e); total += c
        while total > MAX_DEPTH and len(taken) > 1:
            smallest = min(taken, key=est_cov)
            if total - est_cov(smallest) < MIN_DEPTH:
                break
            taken.remove(smallest); total -= est_cov(smallest)
            drop(smallest, "depth cap: dropped smallest cell to bring %.0fx under %g"
                 % (total + est_cov(smallest), MAX_DEPTH))
        if len(taken) < len(cells):
            notes.append("capped at %d of %d cells (~%.0fx)" % (len(taken), len(cells), total))
        note = "; ".join(notes)
        eid = entry_id(sample, chosen)
        for e in taken:
            rows.append([eid, sample, e["prep"], e["name"], "%.1f" % (e["size"] / 1e9),
                         "" if e["cov"] is None else "%.1f" % e["cov"],
                         (note + "; " if note else "") + "prep from " + e["src"]])
        if note:
            notes_out.append("%s: %s" % (eid, note))
    note = " | ".join(notes_out)
    return rows, dropped, note


def refresh(samples, files):
    """Files on S3 not on disk, after the filename-level rules."""
    out = []
    for s in samples:
        on_disk = {movie(e["name"]) for e in files[s]}
        for cohort in COHORTS:
            xml = subprocess.run(["curl", "-s", S3 % (cohort, s)],
                                 capture_output=True, text=True).stdout
            items = re.findall(r"<Key>([^<]+\.bam)</Key><LastModified>[^<]+</LastModified>"
                               r"<ETag>[^<]+</ETag><Size>(\d+)</Size>", xml)
            if not items:
                continue
            new = defaultdict(list)
            for key, size in items:
                name = key.split("/")[-1]
                if "fail_reads" in name or movie(name) in on_disk:
                    continue
                new[movie(name)].append((name, key, int(size)))
            for m, cands in new.items():
                cands.sort(key=lambda c: encoding_rank(c[0]))
                name, key, size = cands[0]
                out.append([s, name, S3_URI % key, size])
            break
    with open(DOWNLOAD, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(out)
    tb = sum(r[3] for r in out) / 1e12
    print("to download: %d files, %.2f TB, %d samples -> %s"
          % (len(out), tb, len({r[0] for r in out}), DOWNLOAD))
    for s in sorted({r[0] for r in out}):
        print("  %-8s %s" % (s, " ".join(r[1] for r in out if r[0] == s)))


def main():
    files = load_files()
    index = {}
    for r in csv.DictReader(open(INDEX)):
        index[r["filename"]] = r
        index.setdefault((r["sample_ID"], movie(r["filename"])), r)
    overrides = load_overrides()
    samples = sorted(files)
    lists = [a for a in sys.argv[1:] if not a.startswith("--")]
    if lists:
        wanted = {l.strip() for l in open(lists[0]) if l.strip() and not l.startswith("#")}
        missing = sorted(wanted - set(samples))
        if missing:
            sys.exit("ERROR: not on disk: %s" % " ".join(missing))
        samples = sorted(wanted)

    plan, excluded, flagged = [], [], []
    for s in samples:
        rows, dropped, note = plan_sample(s, files[s], index, overrides)
        plan += rows
        excluded += dropped
        if note:
            flagged.append((s, note))

    with open(PLAN, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["input_bam", "sample", "prep", "filename", "size_gb", "cov_x", "note"])
        w.writerows(plan)
    with open(EXCLUDED, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "filename", "reason"])
        w.writerows(excluded)

    print("plan: %d files, %d entries, %d samples -> %s"
          % (len(plan), len({r[0] for r in plan}), len({r[1] for r in plan}), PLAN))
    print("excluded: %d files -> %s" % (len(excluded), EXCLUDED))
    reasons = defaultdict(int)
    for _, _, why in excluded:
        reasons[why.split(":")[0].split(" (")[0]] += 1
    for why, n in sorted(reasons.items(), key=lambda x: -x[1]):
        print("  %4d  %s" % (n, why))
    if flagged:
        print("flagged samples:")
        for s, note in flagged:
            print("  %-8s %s" % (s, note))

    if "--refresh" in sys.argv:
        refresh(samples, files)


if __name__ == "__main__":
    main()
