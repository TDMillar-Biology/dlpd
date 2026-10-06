#!/usr/bin/env python3
"""Project reference BED intervals using the existing SVMU synteny mapper.

Requires svmu2 and pandas in the active environment. The BED and delta must use
matching reference sequence names. Input BED is zero-based, half-open; SVMU block
coordinates are assumed to retain NUCMER's one-based base coordinates. We look
up the included endpoint bases start+1 and end, then convert to query BED.

Projection uses block.slope * position + block.y_intercept, exactly as in
dlpd_hic:workflow/scripts/synteny_mapper.py. These are block-based estimates, not gap-aware liftover.
Only uniquely assigned endpoints on the same query sequence with compatible
orientation/order produce a query interval. Missing/ambiguous endpoints remain
in the boundary report. No nearest-block fallback or inward shifting is applied.
"""

from svmu2.orchestration.parse import parse
from svmu2.orchestration.synteny import resolve_synteny
import argparse
from dataclasses import dataclass
from pathlib import Path
import pandas as pd


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delta", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path, help="Reference BED: chrom start end")
    parser.add_argument("--out-prefix", required=True, type=Path)
    args = parser.parse_args()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    return args


def read_intervals(path):
    '''
    read bed format from path as list(interval_objs), ignore comments 
    '''
    intervals = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}:{line} expected chrom start end")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_no}:{line} invalid BED interval")
            intervals.append(Interval(chrom, start, end))
    if not intervals:
        raise ValueError(f"No intervals in {path}")
    return intervals


def boundary_points(intervals):
    ''' 
    list(interval_objs) -> pd.df
    this is for compatibility with dlpd_hic:workflow/scripts/synteny_mapper.py:assign_synteny
    '''
    rows = []
    for interval_id, interval in enumerate(intervals):
        for boundary, position in (("start", interval.start + 1), ("end", interval.end)):
            rows.append(dict(interval_id=interval_id, ref_chr=interval.chrom,
                             ref_start=interval.start, ref_end=interval.end,
                             boundary=boundary, position=position))
    return pd.DataFrame(rows)


def assign_synteny(df, trees):
    """
    dlpd_hic:workflow/scripts/synteny_mapper.py:assign_synteny
    Existing assign_synteny logic, with position replacing midpoint.
    """
    df_out = df.copy()
    block_ids, yhats, query_chroms, strands, statuses = [], [], [], [], []

    for row in df.itertuples(index=False):
        tree = trees.get(row.ref_chr)
        x = row.position
        hits = tree.at(x) if tree is not None else set()

        if len(hits) == 1:
            hit_data = next(iter(hits)).data
            block_ids.append(hit_data.index)
            yhats.append((hit_data.slope * x) + hit_data.y_intercept)
            query_chroms.append(hit_data.query)
            strands.append("+" if hit_data.slope >= 0 else "-")
            statuses.append("mapped")
        else:
            # Match the original mapper: discard missing or nonunique assignments.
            block_ids.append(None)
            yhats.append(None)
            query_chroms.append(None)
            strands.append(None)
            statuses.append("ambiguous" if hits else "unmapped")

    df_out["syntenic_block"] = block_ids
    df_out["yhat"] = yhats
    df_out["query_chr"] = query_chroms
    df_out["strand"] = strands
    df_out["status"] = statuses
    return df_out


def write_outputs(ref_df, out_prefix):
    bed_rows = []
    ref_df = ref_df.copy()
    ref_df["interval_status"] = "unmapped_boundary"
    for _, endpoints in ref_df.groupby("interval_id", sort=False):
        start = endpoints.loc[endpoints["boundary"] == "start"].iloc[0]
        end = endpoints.loc[endpoints["boundary"] == "end"].iloc[0]
        if (endpoints["status"] == "ambiguous").any():
            status = "ambiguous_boundary"
        elif not (endpoints["status"] == "mapped").all():
            status = "unmapped_boundary"
        elif start.query_chr != end.query_chr:
            status = "different_query_sequences"
        elif start.strand != end.strand:
            status = "incompatible_strands"
        else:
            q_start, q_end = int(round(start.yhat)), int(round(end.yhat))
            direction = 1 if start.strand == "+" else -1
            if min(q_start, q_end) < 1 or (
                start.position != end.position and (q_end - q_start) * direction <= 0
            ):
                status = "incompatible_endpoint_order"
            else:
                bed_rows.append((start.query_chr, min(q_start, q_end) - 1, max(q_start, q_end)))
                status = "mapped"
        ref_df.loc[endpoints.index, "interval_status"] = status

    with Path(str(out_prefix) + ".euchromatin.bed").open("w") as handle:
        for chrom, start, end in bed_rows:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    ref_df.to_csv(str(out_prefix) + ".boundaries.tsv", sep="\t", index=False, na_rep=".")
    print(f"Mapped {len(bed_rows)}/{ref_df['interval_id'].nunique()} intervals")


def main():
    args = parse_args()
    REF_euchromatin = read_intervals(args.bed)

    _, primary_alns = parse(args.delta)
    resolve_synteny(primary_alignments=primary_alns, breakpoint_map=None)

    ref_trees = {aln.reference: aln.reference_synteny_tree for aln in primary_alns.values()}

    ref_df = assign_synteny(boundary_points(REF_euchromatin), ref_trees)
    write_outputs(ref_df, args.out_prefix)


if __name__ == "__main__":
    main()
