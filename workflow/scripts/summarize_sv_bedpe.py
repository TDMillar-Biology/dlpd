#!/usr/bin/env python3
"""
Summarize SVMU2 BEDPE structural-variant calls across samples.

Expected modern SVMU2 BEDPE format (11 columns):
    ref_chrom ref_start ref_end query_chrom query_start query_end
    event_id score ref_strand query_strand sv_type

Legacy BEDPE6 files are accepted, but sv_type is recorded as UNKNOWN and
they are excluded from SV-type plots unless --include-unknown is used.

Outputs:
    events.tsv
    sample_svtype_summary.tsv
    sample_chrom_svtype_summary.tsv
    sample_summary.tsv
    sv_counts_by_sample.png
    ref_bp_by_sample.png
    query_bp_by_sample.png
    delta_bp_by_sample.png
    chromosome_ref_bp.png
    chromosome_query_bp.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
import matplotlib.pyplot as plt


BEDPE6 = [
    "ref_chrom", "ref_start", "ref_end",
    "query_chrom", "query_start", "query_end",
]

BEDPE11 = BEDPE6 + [
    "event_id", "score", "ref_strand", "query_strand", "sv_type",
]

SV_ORDER = ["INS", "DEL", "DUP", "DUP_TANDEM", "INV", "COMPLEX", "UNKNOWN"]
CHROM_ORDER = ["X", "2L", "2R", "3L", "3R", "4", "Y", "mitochondrion_genome"]


def sample_from_path(path: Path) -> str:
    """Derive sample name from '<sample>.svmu2.bedpe'."""
    name = path.name
    suffix = ".svmu2.bedpe"
    return name[:-len(suffix)] if name.endswith(suffix) else path.stem


def read_bedpe(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")

    if df.shape[1] == 11:
        df.columns = BEDPE11
        modern = True
    elif df.shape[1] == 6:
        df.columns = BEDPE6
        df["event_id"] = [f"{sample_from_path(path)}_{i:06d}" for i in range(len(df))]
        df["score"] = "."
        df["ref_strand"] = "."
        df["query_strand"] = "."
        df["sv_type"] = "UNKNOWN"
        modern = False
    else:
        raise ValueError(
            f"{path}: expected 6 or 11 tab-separated columns, found {df.shape[1]}"
        )

    for col in ["ref_start", "ref_end", "query_start", "query_end"]:
        df[col] = pd.to_numeric(df[col], errors="raise").astype("int64")

    # Defensive normalization.
    r0 = df[["ref_start", "ref_end"]].min(axis=1)
    r1 = df[["ref_start", "ref_end"]].max(axis=1)
    q0 = df[["query_start", "query_end"]].min(axis=1)
    q1 = df[["query_start", "query_end"]].max(axis=1)

    df["ref_start"] = r0
    df["ref_end"] = r1
    df["query_start"] = q0
    df["query_end"] = q1

    df["ref_span_bp"] = df["ref_end"] - df["ref_start"]
    df["query_span_bp"] = df["query_end"] - df["query_start"]

    # Query allele span minus reference allele span.
    # Useful as a signed length-change descriptor for INS/DEL.
    # Do not interpret this automatically as copy-number change for
    # DUP/DUP_TANDEM/INV/COMPLEX events.
    df["span_delta_bp"] = df["query_span_bp"] - df["ref_span_bp"]
    df["abs_span_delta_bp"] = df["span_delta_bp"].abs()

    df["sample"] = sample_from_path(path)
    df["source_file"] = str(path)
    df["bedpe_format"] = "BEDPE11" if modern else "BEDPE6"

    return df


def discover_files(inputs: list[str], root: str | None) -> list[Path]:
    paths = [Path(x) for x in inputs]

    if root:
        paths.extend(Path(root).rglob("*.svmu2.bedpe"))

    files = sorted({p.resolve() for p in paths if p.is_file()})

    if not files:
        raise FileNotFoundError(
            "No BEDPE files found. Supply files directly or use --root."
        )

    return files


def ordered_columns_present(values, preferred):
    values = list(pd.unique(values))
    return [x for x in preferred if x in values] + sorted(
        x for x in values if x not in preferred
    )


def make_stacked_sample_plot(
    summary: pd.DataFrame,
    value_col: str,
    ylabel: str,
    output: Path,
):
    pivot = summary.pivot_table(
        index="sample",
        columns="sv_type",
        values=value_col,
        aggfunc="sum",
        fill_value=0,
    )

    order = [x for x in SV_ORDER if x in pivot.columns]
    pivot = pivot.reindex(columns=order)

    ax = pivot.plot(kind="bar", stacked=True, figsize=(11, 6))
    ax.set_xlabel("Sample")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " by sample and SV type")
    ax.legend(title="SV type", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def make_chromosome_plot(
    summary: pd.DataFrame,
    value_col: str,
    ylabel: str,
    output: Path,
):
    # Sum across samples to show the chromosome-level composition of burden.
    grouped = (
        summary.groupby(["ref_chrom", "sv_type"], as_index=False)[value_col]
        .sum()
    )

    pivot = grouped.pivot_table(
        index="ref_chrom",
        columns="sv_type",
        values=value_col,
        aggfunc="sum",
        fill_value=0,
    )

    chroms = ordered_columns_present(pivot.index, CHROM_ORDER)
    types = [x for x in SV_ORDER if x in pivot.columns]
    pivot = pivot.reindex(index=chroms, columns=types)

    ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6))
    ax.set_xlabel("Reference chromosome")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " by chromosome and SV type")
    ax.legend(title="SV type", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bedpe",
        nargs="*",
        help="BEDPE files to process.",
    )
    parser.add_argument(
        "--root",
        help="Recursively search this directory for *.svmu2.bedpe.",
    )
    parser.add_argument(
        "--outdir",
        default="sv_burden_summary",
        help="Output directory (default: sv_burden_summary).",
    )
    parser.add_argument(
        "--include-unknown",
        action="store_true",
        help="Include legacy BEDPE6 UNKNOWN calls in SV-type plots.",
    )
    args = parser.parse_args()

    files = discover_files(args.bedpe, args.root)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    frames = []
    for path in files:
        df = read_bedpe(path)
        frames.append(df)
        print(
            f"{df['sample'].iat[0]}: {len(df):,} events "
            f"({df['bedpe_format'].iat[0]})"
        )

    events = pd.concat(frames, ignore_index=True)

    # Per-event table: this is the canonical tidy dataset for plotting/reanalysis.
    event_cols = [
        "sample", "ref_chrom", "ref_start", "ref_end",
        "query_chrom", "query_start", "query_end",
        "event_id", "ref_strand", "query_strand", "sv_type",
        "ref_span_bp", "query_span_bp", "span_delta_bp",
        "abs_span_delta_bp", "bedpe_format", "source_file",
    ]
    events[event_cols].to_csv(
        outdir / "events.tsv", sep="\t", index=False
    )

    sample_svtype = (
        events.groupby(["sample", "sv_type"], as_index=False)
        .agg(
            n_events=("event_id", "size"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            net_span_delta_bp=("span_delta_bp", "sum"),
            abs_span_delta_bp=("abs_span_delta_bp", "sum"),
            median_ref_span_bp=("ref_span_bp", "median"),
            median_query_span_bp=("query_span_bp", "median"),
        )
    )
    sample_svtype.to_csv(
        outdir / "sample_svtype_summary.tsv", sep="\t", index=False
    )

    sample_chrom_svtype = (
        events.groupby(["sample", "ref_chrom", "sv_type"], as_index=False)
        .agg(
            n_events=("event_id", "size"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            net_span_delta_bp=("span_delta_bp", "sum"),
            abs_span_delta_bp=("abs_span_delta_bp", "sum"),
        )
    )
    sample_chrom_svtype.to_csv(
        outdir / "sample_chrom_svtype_summary.tsv", sep="\t", index=False
    )

    sample_summary = (
        events.groupby("sample", as_index=False)
        .agg(
            n_events=("event_id", "size"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            net_span_delta_bp=("span_delta_bp", "sum"),
            abs_span_delta_bp=("abs_span_delta_bp", "sum"),
        )
    )
    sample_summary.to_csv(
        outdir / "sample_summary.tsv", sep="\t", index=False
    )

    plot_events = events
    if not args.include_unknown:
        plot_events = plot_events[plot_events["sv_type"] != "UNKNOWN"]

    if plot_events.empty:
        print(
            "\nNo typed BEDPE11 events are available, so plots were not made.\n"
            "Regenerate BEDPE files with the 11-column exporter, then rerun."
        )
        return

    plot_sample_svtype = (
        plot_events.groupby(["sample", "sv_type"], as_index=False)
        .agg(
            n_events=("event_id", "size"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
            net_span_delta_bp=("span_delta_bp", "sum"),
        )
    )

    plot_chrom = (
        plot_events.groupby(["sample", "ref_chrom", "sv_type"], as_index=False)
        .agg(
            n_events=("event_id", "size"),
            ref_bp_affected=("ref_span_bp", "sum"),
            query_bp_affected=("query_span_bp", "sum"),
        )
    )

    make_stacked_sample_plot(
        plot_sample_svtype,
        "n_events",
        "SV count",
        outdir / "sv_counts_by_sample.png",
    )
    make_stacked_sample_plot(
        plot_sample_svtype,
        "ref_bp_affected",
        "Reference-space bp affected",
        outdir / "ref_bp_by_sample.png",
    )
    make_stacked_sample_plot(
        plot_sample_svtype,
        "query_bp_affected",
        "Query-space bp affected",
        outdir / "query_bp_by_sample.png",
    )
    make_stacked_sample_plot(
        plot_sample_svtype,
        "net_span_delta_bp",
        "Net query-minus-reference span (bp)",
        outdir / "delta_bp_by_sample.png",
    )
    make_chromosome_plot(
        plot_chrom,
        "ref_bp_affected",
        "Reference-space bp affected",
        outdir / "chromosome_ref_bp.png",
    )
    make_chromosome_plot(
        plot_chrom,
        "query_bp_affected",
        "Query-space bp affected",
        outdir / "chromosome_query_bp.png",
    )

    print(f"\nWrote summaries and plots to: {outdir}")


if __name__ == "__main__":
    main()
