#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_geometry(path):
    return pd.read_csv(path, sep="\t")


def analyze_quotient_thresholds(df, thresholds):
    complex_df = df[df["sv_type"] == "COMPLEX"].copy()

    rows = []

    for threshold in thresholds:
        deletion_like = (
            (complex_df["ref_span_bp"] > complex_df["query_span_bp"])
            & (
                complex_df["query_span_bp"]
                / complex_df["ref_span_bp"]
                <= threshold
            )
        )

        insertion_like = (
            (complex_df["query_span_bp"] > complex_df["ref_span_bp"])
            & (
                complex_df["ref_span_bp"]
                / complex_df["query_span_bp"]
                <= threshold
            )
        )

        n_complex = len(complex_df)
        n_del = int(deletion_like.sum())
        n_ins = int(insertion_like.sum())
        n_simple_like = n_del + n_ins

        rows.append(
            {
                "sample": complex_df["sample"].iloc[0],
                "threshold": threshold,
                "n_complex": n_complex,
                "n_deletion_like": n_del,
                "n_insertion_like": n_ins,
                "n_simple_like": n_simple_like,
                "fraction_simple_like": (
                    n_simple_like / n_complex
                    if n_complex
                    else 0
                ),
            }
        )

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Quotient-threshold analysis of COMPLEX SVMU2 calls."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="*.sv_geometry.tsv from summarize_svs.py.",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV.",
    )

    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=[0.01, 0.05, 0.10],
        help="Minor/major span-ratio thresholds.",
    )

    args = parser.parse_args()

    df = parse_geometry(args.input)
    summary = analyze_quotient_thresholds(df, args.thresholds)
    summary.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
