#!/usr/bin/env python3

'''
Example usage: 
python downsample_basesSequenced.py \
  --root . \
  --barcodes 1-6

Or to enforce a target
python downsample_bams_by_bases.py \
  --root . \
  --barcodes 1-6 \ 
  --target-bases NUMBER 
'''
#!/usr/bin/env python3

import argparse
import random
import sys
from pathlib import Path

import pysam


def parse_barcode_selection(selection):
    """
    Parse barcode selection strings like:
      "1-6"
      "1,2,5,7"
      "1-3,8,10-12"

    Returns a set of integers, or None if no selection was provided.
    """
    if selection is None:
        return None

    selected = set()

    for part in selection.split(","):
        part = part.strip()
        if not part:
            continue

        if "-" in part:
            start, end = part.split("-", 1)
            start = int(start)
            end = int(end)
            if start > end:
                raise ValueError(f"Invalid barcode range: {part}")
            for i in range(start, end + 1):
                selected.add(i)
        else:
            selected.add(int(part))

    return selected


def get_read_length(read):
    """
    Return read/query length robustly for aligned or unaligned BAMs.
    """
    if read.query_sequence is not None:
        return len(read.query_sequence)

    if read.query_length is not None:
        return read.query_length

    return 0


def find_input_bam(barcode_dir):
    """
    Find the unaligned BAM in a barcode directory.

    Rule:
      - Must end with .bam
      - Must start with 'PAQ'

    Returns:
      Path to BAM, or None if not found
    """
    candidate_bams = sorted(
        [
            p for p in barcode_dir.iterdir()
            if p.is_file() and p.suffix == ".bam" and p.name.startswith("PAQ")
        ]
    )

    if len(candidate_bams) == 0:
        return None

    if len(candidate_bams) > 1:
        print(
            f"[WARN] Multiple PAQ BAMs found in {barcode_dir}; using first: {candidate_bams[0].name}",
            file=sys.stderr
        )

    return candidate_bams[0]


def find_barcode_bams(root_dir, selected_barcodes=None):
    """
    Structure like: 
      root/
        NoAdaptive/
          demux/
            barcode01/
              PAQ....bam
            barcode02/
              PAQ....bam
        Pmand_enriched/
          demux/
            barcode01/
              PAQ....bam
        Stub_depleted/
          demux/
            barcode01/
              PAQ....bam

    Returns:
      list of dicts with keys:
        type_name, barcode_name, barcode_num, bam_path, barcode_dir
    """
    root = Path(root_dir).resolve()
    results = []

    for type_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        demux_dir = type_dir / "demux"
        if not demux_dir.is_dir():
            continue

        barcode_dirs = sorted(
            [p for p in demux_dir.iterdir() if p.is_dir() and p.name.startswith("barcode")]
        )

        for barcode_dir in barcode_dirs:
            try:
                barcode_num = int(barcode_dir.name.replace("barcode", ""))
            except ValueError:
                print(f"[WARN] Skipping unexpected barcode directory name: {barcode_dir}", file=sys.stderr)
                continue

            if selected_barcodes is not None and barcode_num not in selected_barcodes:
                continue

            bam_path = find_input_bam(barcode_dir)

            if bam_path is None:
                print(
                    f"[WARN] No PAQ*.bam unaligned BAM found in {barcode_dir}",
                    file=sys.stderr
                )
                continue

            results.append(
                {
                    "type_name": type_dir.name,
                    "barcode_name": barcode_dir.name,
                    "barcode_num": barcode_num,
                    "bam_path": bam_path,
                    "barcode_dir": barcode_dir,
                }
            )

    return results


def collect_bam_stats(bam_path):
    """
    Return per-read lengths plus summary stats.
    """
    read_lengths = []
    total_reads = 0
    total_bases = 0

    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for read in bam:
            read_len = get_read_length(read)
            read_lengths.append(read_len)
            total_reads += 1
            total_bases += read_len

    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "read_lengths": read_lengths,
    }


def choose_reads_by_base_target(read_lengths, target_bases, seed=1):
    """
    Downsample by:

    - randomize read order
    - only add the next read if adding it makes total bases closer to target

    Returns:
      selected_indices (set of original read indices),
      selected_reads_count,
      selected_bases
    """
    rng = random.Random(seed)

    indexed_lengths = list(enumerate(read_lengths))
    rng.shuffle(indexed_lengths)

    selected_indices = set()
    current_bases = 0
    current_reads = 0

    for idx, read_len in indexed_lengths:
        before = abs(target_bases - current_bases)
        after = abs(target_bases - (current_bases + read_len))

        if after < before:
            selected_indices.add(idx)
            current_bases += read_len
            current_reads += 1

            if current_bases == target_bases:
                break

    return selected_indices, current_reads, current_bases


def write_downsampled_bam(input_bam, output_bam, selected_indices):
    """
    Write selected reads by original read index order.
    """
    written = 0

    with pysam.AlignmentFile(str(input_bam), "rb", check_sq=False) as in_bam:
        with pysam.AlignmentFile(str(output_bam), "wb", template=in_bam) as out_bam:
            for idx, read in enumerate(in_bam):
                if idx in selected_indices:
                    out_bam.write(read)
                    written += 1

    return written


def format_int(x):
    return f"{x:,}"


def print_table(rows, title=None):
    if title:
        print(f"\n{title}")

    headers = [
        "type",
        "barcode",
        "input_bam",
        "orig_reads",
        "orig_bases",
        "down_reads",
        "down_bases",
        "diff_from_target",
        "output_bam",
    ]

    widths = {h: len(h) for h in headers}
    for row in rows:
        for h in headers:
            widths[h] = max(widths[h], len(str(row.get(h, ""))))

    header_line = "  ".join(h.ljust(widths[h]) for h in headers)
    sep_line = "  ".join("-" * widths[h] for h in headers)

    print(header_line)
    print(sep_line)

    for row in rows:
        print("  ".join(str(row.get(h, "")).ljust(widths[h]) for h in headers))


def save_stats_tsv(rows, out_path):
    headers = [
        "type",
        "barcode",
        "input_bam",
        "orig_reads",
        "orig_bases",
        "down_reads",
        "down_bases",
        "diff_from_target",
        "output_bam",
    ]

    with open(out_path, "w") as out:
        out.write("\t".join(headers) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Downsample unaligned PAQ*.bam files across barcode folders to a shared base-count target.\n"
            "The shared target is the minimum total bases across all selected BAMs unless overridden."
        )
    )
    parser.add_argument(
        "--root",
        required=True,
        help="Root directory containing NoAdaptive/, Pmand_enriched/, Stub_depleted/, etc."
    )
    parser.add_argument(
        "--barcodes",
        help="Barcode selection, e.g. 1-6 or 1,2,5,7. Default: all barcodes."
    )
    parser.add_argument(
        "--target-bases",
        type=int,
        default=None,
        help="Override shared target base count. Default: minimum total bases across selected files."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for reproducible downsampling. Default: 1"
    )
    parser.add_argument(
        "--suffix",
        default="downsampled",
        help="Suffix to use in output BAM names. Default: downsampled"
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite downsampled BAMs if they already exist."
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    if not root.is_dir():
        print(f"[ERROR] Root directory does not exist: {root}", file=sys.stderr)
        sys.exit(1)

    try:
        selected_barcodes = parse_barcode_selection(args.barcodes)
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    jobs = find_barcode_bams(root, selected_barcodes=selected_barcodes)
    if not jobs:
        print("[ERROR] No matching PAQ*.bam files found.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Found {len(jobs)} PAQ BAM file(s).")
    print("[INFO] Collecting original read/base statistics...")

    for job in jobs:
        stats = collect_bam_stats(job["bam_path"])
        job["orig_reads"] = stats["total_reads"]
        job["orig_bases"] = stats["total_bases"]
        job["read_lengths"] = stats["read_lengths"]

    if args.target_bases is not None:
        shared_target = args.target_bases
        print(f"\n[INFO] Using user-specified shared target: {format_int(shared_target)} bases")
    else:
        shared_target = min(job["orig_bases"] for job in jobs)
        print(f"\n[INFO] Computed shared target (minimum total bases): {format_int(shared_target)} bases")

    pre_rows = []
    for job in jobs:
        pre_rows.append(
            {
                "type": job["type_name"],
                "barcode": job["barcode_name"],
                "input_bam": job["bam_path"].name,
                "orig_reads": format_int(job["orig_reads"]),
                "orig_bases": format_int(job["orig_bases"]),
                "down_reads": "",
                "down_bases": "",
                "diff_from_target": "",
                "output_bam": "",
            }
        )

    print_table(pre_rows, title="Original BAM statistics")

    print("\n[INFO] Downsampling BAMs...")

    final_rows = []

    for job in jobs:
        type_name = job["type_name"]
        barcode_name = job["barcode_name"]
        bam_path = job["bam_path"]
        barcode_dir = job["barcode_dir"]

        output_bam = barcode_dir / f"{type_name}-{barcode_name}.{args.suffix}.{shared_target}bp.bam"

        if output_bam.exists() and not args.overwrite:
            print(f"[SKIP] Output exists: {output_bam}")
            out_stats = collect_bam_stats(output_bam)
            down_reads = out_stats["total_reads"]
            down_bases = out_stats["total_bases"]
        else:
            selected_indices, down_reads, down_bases = choose_reads_by_base_target(
                job["read_lengths"],
                shared_target,
                seed=args.seed
            )

            print(
                f"[RUN] {bam_path.name} -> {output_bam.name} | "
                f"orig={format_int(job['orig_bases'])} bp, "
                f"target={format_int(shared_target)} bp, "
                f"down={format_int(down_bases)} bp"
            )

            written = write_downsampled_bam(bam_path, output_bam, selected_indices)

            if written != down_reads:
                print(
                    f"[WARN] Expected to write {down_reads} reads but wrote {written} reads for {bam_path.name}",
                    file=sys.stderr
                )
                down_reads = written

        diff_from_target = abs(shared_target - down_bases)

        final_rows.append(
            {
                "type": type_name,
                "barcode": barcode_name,
                "input_bam": bam_path.name,
                "orig_reads": format_int(job["orig_reads"]),
                "orig_bases": format_int(job["orig_bases"]),
                "down_reads": format_int(down_reads),
                "down_bases": format_int(down_bases),
                "diff_from_target": format_int(diff_from_target),
                "output_bam": output_bam.name,
            }
        )

    print_table(final_rows, title="Downsampled BAM statistics")

    stats_tsv = root / f"downsampling_stats.{shared_target}bp.tsv"
    save_stats_tsv(final_rows, stats_tsv)
    print(f"\n[INFO] Saved stats table to: {stats_tsv}")
    print(f"[INFO] Shared target used: {format_int(shared_target)} bases")


if __name__ == "__main__":
    main()