#!/usr/bin/env python3

import argparse

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


SV_ORDER = [
    "INS",
    "DEL",
    "DUP",
    "DUP_TANDEM",
    "INV",
    "COMPLEX",
]


THETA_BOUNDARIES = [-135, -90, -45, 45, 90, 135]


def parse_geometry(path):
    return pd.read_csv(path, sep="\t")


def filter_sv_types(df, ignored_types):
    if not ignored_types:
        return df

    return df[~df["sv_type"].isin(ignored_types)].copy()


def plot_span_geometry(df, output):
    plt.figure(figsize=(7, 7))

    for sv_type in SV_ORDER:
        subset = df[df["sv_type"] == sv_type]

        if subset.empty:
            continue

        plt.scatter(
            subset["ref_span_bp"],
            subset["query_span_bp"],
            s=14,
            alpha=0.5,
            label=sv_type,
        )

    max_span = max(
        df["ref_span_bp"].max(),
        df["query_span_bp"].max(),
    )

    plt.plot(
        [1, max_span],
        [1, max_span],
        linestyle="--",
        linewidth=1,
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Reference span (bp)")
    plt.ylabel("Query span (bp)")
    plt.title(f"{df['sample'].iloc[0]} SV geometry")
    plt.legend(title="SV type")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def plot_theta_distribution(df, output):
    plt.figure(figsize=(10, 5))

    for sv_type in SV_ORDER:
        subset = df[df["sv_type"] == sv_type]

        if subset.empty:
            continue

        plt.hist(
            subset["theta"],
            bins=36,
            alpha=0.5,
            label=sv_type,
        )

    for theta in THETA_BOUNDARIES:
        plt.axvline(
            theta,
            linestyle="--",
            linewidth=1,
        )

    plt.xticks(
        [-180, -135, -90, -45, 0, 45, 90, 135, 180]
    )

    plt.xlabel("Segment angle (degrees)")
    plt.ylabel("Number of SVs")
    plt.title(f"{df['sample'].iloc[0]} SV angle distribution")
    plt.legend(title="SV type")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot SVMU2 structural variant geometry by SV type."
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
        "--ignore-sv-type",
        action="append",
        default=[],
        help=(
            "SV type to exclude from plots. "
            "May be supplied more than once."
        ),
    )

    args = parser.parse_args()

    df = parse_geometry(args.input)
    df = filter_sv_types(df, args.ignore_sv_type)

    if df.empty:
        raise ValueError("No SVs remain after filtering.")

    plot_span_geometry(
        df,
        f"{args.out_prefix}.span_geometry.png",
    )

    plot_theta_distribution(
        df,
        f"{args.out_prefix}.theta_distribution.png",
    )


if __name__ == "__main__":
    main()
