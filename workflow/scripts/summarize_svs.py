#!/usr/bin/env python3

import argparse
import math
import pandas as pd

BEDPE11_COLUMNS = [
    "ref_chrom", "ref_start", "ref_end",
    "query_chrom", "query_start", "query_end",
    "event_id", "score", "ref_strand", "query_strand", "sv_type",
]

BEDPE12_COLUMNS = BEDPE11_COLUMNS + ["theta"]


def parse_bedpe(path, sample):
    df = pd.read_csv(path, sep="\t", header=None)

    if df.shape[1] == 11:
        df.columns = BEDPE11_COLUMNS
    elif df.shape[1] == 12:
        df.columns = BEDPE12_COLUMNS
    else:
        raise ValueError(
            f"Expected BEDPE11 or BEDPE12 input, found {df.shape[1]} columns."
        )

    df["sample"] = sample
    return df


def add_geometry_columns(df):
    df = df.copy()

    df["ref_span_bp"] = df["ref_end"] - df["ref_start"]
    df["query_span_bp"] = df["query_end"] - df["query_start"]
    df["span_delta_bp"] = df["query_span_bp"] - df["ref_span_bp"]
    df["abs_span_delta_bp"] = df["span_delta_bp"].abs()

    max_span = df[["ref_span_bp", "query_span_bp"]].max(axis=1)
    min_span = df[["ref_span_bp", "query_span_bp"]].min(axis=1)

    df["span_ratio"] = (
        df["query_span_bp"] / df["ref_span_bp"].replace(0, pd.NA)
    )
    df["span_similarity"] = min_span / max_span.replace(0, pd.NA)

    if "theta" not in df.columns:
        ref_direction = df["ref_strand"].map({"+": 1, "-": -1})
        query_direction = df["query_strand"].map({"+": 1, "-": -1})

        dx = df["ref_span_bp"] * ref_direction
        dy = df["query_span_bp"] * query_direction

        df["theta"] = [
            round(math.degrees(math.atan2(y, x)), 1)
            for x, y in zip(dx, dy)
        ]
    else:
        df["theta"] = df["theta"].round(1)

    return df


def build_geometry_table(df):
    return df[[
        "sample",
        "ref_chrom",
        "ref_start",
        "ref_end",
        "query_chrom",
        "query_start",
        "query_end",
        "event_id",
        "theta",
        "score",
        "ref_strand",
        "query_strand",
        "sv_type",
        "ref_span_bp",
        "query_span_bp",
        "span_delta_bp",
        "abs_span_delta_bp",
        "span_ratio",
        "span_similarity",
    ]].copy()


def summarize_by_chromosome(df):
    return (
        df.groupby(["sample", "ref_chrom", "sv_type"])
        .agg(
            n_events=("event_id", "count"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            span_delta_bp=("span_delta_bp", "sum"),
            median_ref_span_bp=("ref_span_bp", "median"),
            median_query_span_bp=("query_span_bp", "median"),
        )
        .reset_index()
    )


def summarize_by_sample(df):
    return (
        df.groupby(["sample", "sv_type"])
        .agg(
            n_events=("event_id", "count"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            span_delta_bp=("span_delta_bp", "sum"),
            median_ref_span_bp=("ref_span_bp", "median"),
            median_query_span_bp=("query_span_bp", "median"),
        )
        .reset_index()
    )


def write_tsv(df, path):
    df.to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Summarize SVMU2 BEDPE structural variant calls."
    )
    parser.add_argument("--bedpe", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    bedpe = parse_bedpe(args.bedpe, args.sample)
    bedpe = add_geometry_columns(bedpe)

    geometry = build_geometry_table(bedpe)
    chromosome_summary = summarize_by_chromosome(bedpe)
    sample_summary = summarize_by_sample(bedpe)

    write_tsv(geometry, f"{args.out_prefix}.sv_geometry.tsv")
    write_tsv(
        chromosome_summary,
        f"{args.out_prefix}.chromosome_summary.tsv",
    )
    write_tsv(
        sample_summary,
        f"{args.out_prefix}.sample_summary.tsv",
    )


if __name__ == "__main__":
    main()
