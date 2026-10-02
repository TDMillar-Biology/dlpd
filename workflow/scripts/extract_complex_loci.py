#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_geometry(path):
    return pd.read_csv(path, sep="\t")


def make_bed(df, chrom_col, start_col, end_col, strand_col, padding):
    bed = pd.DataFrame()

    bed["chrom"] = df[chrom_col]
    bed["start"] = (df[start_col] - padding).clip(lower=0)
    bed["end"] = df[end_col] + padding
    bed["name"] = df["event_id"]
    bed["score"] = 0
    bed["strand"] = df[strand_col]

    return bed


def main():
    parser = argparse.ArgumentParser(
        description="Create reference/query BED files for COMPLEX SVMU2 calls."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="*.sv_geometry.tsv from summarize_svs.py.",
    )

    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix.",
    )

    parser.add_argument(
        "--padding",
        type=int,
        default=500,
        help="Flanking sequence to add on each side (default: 500 bp).",
    )

    args = parser.parse_args()

    df = parse_geometry(args.input)
    complex_df = df[df["sv_type"] == "COMPLEX"].copy()

    ref_bed = make_bed(
        complex_df,
        "ref_chrom",
        "ref_start",
        "ref_end",
        "ref_strand",
        args.padding,
    )

    query_bed = make_bed(
        complex_df,
        "query_chrom",
        "query_start",
        "query_end",
        "query_strand",
        args.padding,
    )

    ref_bed.to_csv(
        f"{args.out_prefix}.complex.reference.bed",
        sep="\t",
        header=False,
        index=False,
    )

    query_bed.to_csv(
        f"{args.out_prefix}.complex.query.bed",
        sep="\t",
        header=False,
        index=False,
    )


if __name__ == "__main__":
    main()
