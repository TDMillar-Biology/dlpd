"""
delta_parser.py

Current as of: Jan 15, 2025

Trevor D. Millar

Core utilities for parsing MUMmer delta files and defining
alignment blocks and primary synteny intervals.

Synteny is defined at the alignment-block level and is not
assumed to be base-pair exact.
"""

import sys
import os
import argparse
from matplotlib import pyplot as plt
from more_itertools import peekable
from collections import defaultdict
from intervaltree import IntervalTree
from dataclasses import dataclass
from math import sqrt

@dataclass(frozen=True)
class BoundingBox:
    ''' 
    Holds the most extreme coordinates (x,y) of rectangle containg the minimal 
    space required to hold every alignment block in a 2D plane. Is immutable
    '''
    bottom_left: tuple
    top_left: tuple
    top_right: tuple
    bottom_right: tuple

class Alignment:
    ''' 
    This is the full alignment reported by the aligner for a single chromosome
    It is composed of multiple individual line segments denoting synteny
    '''
    def __init__(self, reference, query, reference_length, query_length):
        self.reference = reference
        self.query = query
        self.reference_length = int(reference_length)
        self.query_length = int(query_length)
        self.alignment_blocks = None
        self.zero_coverage_sites = 0
        self.covered_sites = 0
        self.coverage_by_position = []
        self.traversal_len = None
        self.origin = (0,0) ## poised for removal
        self.slope = None ## poised for removal
        self.primary_synteny_blocks = None # Djikstra traversal defines these
        self.primary_synteny_tree = None # main diagonal / orthologous
        self.exhaustive_tree = None   # all homologous blocks

    def __str__(self):
        return f'{self.reference} by {self.query} alignment, consisting of {0 if self.alignment_blocks is None else len(self.alignment_blocks)} alignment blocks'

    def add_alignment_block(self, line, indel_map, index):
        if self.alignment_blocks is None:
            self.alignment_blocks = []
        self.alignment_blocks.append(AlignmentBlock(self.reference, *line, self.query, indel_map, index))
        self.alignment_blocks.sort(key=lambda x: x.left_most) ## keep the alignment blocks sorted by leftmost position on reference

    def set_primary_synteny_blocks(self, blocks):
        """
        Register primary syntenic alignment blocks
        (e.g. result of Dijkstra traversal).
        """
        if blocks is None or len(blocks) == 0:
            raise ValueError("Primary synteny block list is empty")

        self.primary_synteny_blocks = blocks
        self.primary_synteny_tree = None  # invalidate any stale tree

    def calculate_coverage(self, ref_len):
        ''' 
        Treat start,end of alignment as events with value (+1,-1)
        Find segments along ref length that have no alignments
        '''
        events = []
        for block in self.alignment_blocks:
            events.append((block.left_most, 1))
            events.append((block.right_most + 1, -1))
        
        events.sort()
        coverage = 0
        last_pos = 0
        self.zero_coverage_sites = 0
        
        for pos, change in events: # go through events, add to zero_coverage_sites the uncovered region
            if pos > last_pos and coverage == 0:
                self.zero_coverage_sites += pos - last_pos  # count uncovered region
            coverage += change
            last_pos = pos

        # Count trailing uncovered region
        if last_pos < ref_len:
            self.zero_coverage_sites += ref_len - last_pos

        self.covered_sites = ref_len - self.zero_coverage_sites

    def build_exhaustive_tree(self):
        from intervaltree import IntervalTree

        if self.alignment_blocks is None:
            raise RuntimeError("No alignment blocks loaded")

        tree = IntervalTree()
        for block in self.alignment_blocks:
            tree[block.left_most:block.right_most] = block

        self.exhaustive_tree = tree

    def build_primary_synteny_tree(self):
        if self.primary_synteny_blocks is None:
            raise RuntimeError("primary synteny blocks not set, use djikstra traversal before this tool")
        if len(self.primary_synteny_blocks) == 0:
            raise RuntimeError("primary synteny block list is empty (this is trouble)")
        from intervaltree import IntervalTree
        tree = IntervalTree()
        for block in self.primary_synteny_blocks:
            tree[block.left_most:block.right_most] = block
        self.primary_synteny_tree = tree

    def build_synteny(self, max_jump=100_000, weights="euclidean"):
        """
        Fully prepare alignment for synteny queries.
        """

        # 1) Exhaustive homologous blocks
        if self.exhaustive_tree is None:
            self.build_exhaustive_tree()

        # 2) Trend + geometry
        self.evaluate_alignment_trend()
        self.compute_bounding_box()
        self.find_source_sink_nodes()

        # 3) Primary synteny traversal
        from svmu_graph import compute_primary_synteny
        compute_primary_synteny(self, max_jump=max_jump, weights=weights)


    def compute_bounding_box(self):
        """
        Computes an immutable BoundingBox to the alignment object (see BoundingBox class definition)
        """
        if not self.alignment_blocks:
            raise ValueError("No alignment blocks provided.")

        minx = min(b.reference_start for b in self.alignment_blocks)
        maxx = max(b.reference_end for b in self.alignment_blocks)
        miny = min(b.query_start for b in self.alignment_blocks)
        maxy = max(b.query_end for b in self.alignment_blocks)

        bottom_left  = (minx, miny)
        top_left     = (minx, maxy)
        top_right    = (maxx, maxy)
        bottom_right = (maxx, miny)

        self.bounding_box = BoundingBox(
            bottom_left=bottom_left,
            top_left=top_left,
            top_right=top_right,
            bottom_right=bottom_right
        )
    @staticmethod
    def euclidean_distance(p1, p2):
        '''Need to be careful here about namespace, would be good to check, maybe define globally if possible'''
        if len(p1) != len(p2):
            raise ValueError("Points must have the same dimensionality")
        return sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
    
    @staticmethod
    def x_distance(p1, p2):
        """
        Reference-space (x-axis) distance only.
        Used when reference order should dominate (primary synteny endpoints).
        """
        if len(p1) < 1 or len(p2) < 1:
            raise ValueError("Points must have at least one dimension")

        return abs(p2[0] - p1[0])

    def evaluate_alignment_trend(self, debug_dir="svmu_debugging"):
        """
        Analyze monotonic trend in cumulative block lengths.
        
        source and sink node discovery has been moved
        to function find_source_sink_nodes()

        Returns
        None
            Alignment is monotonic (expected case)
        dict
            Diagnostic information for non-monotonic alignments
        """
        import numpy as np
        import os
        import matplotlib.pyplot as plt
        import pymannkendall as mk

        # trace of the weighted length of the alignment (weighted length = slope * length)
        trace = [b.weighted_length for b in self.alignment_blocks]
        cumsum = np.cumsum(trace)

        if len(cumsum) < 3:
            # Too few points to assess trend meaningfully
            # We may need to throw an exception here... undetermined
            return None
        
        # Read more on this, determines if a function is monotonic decrease / increase
        # works well for a quick check for monotonicity -- may need to delve deeper into the stats here
        result = mk.hamed_rao_modification_test(cumsum) 
        
        self.trend = result.trend # result.trend belongs to {"increasing", "decreasing", "no trend"}
        self.trend_checked = True
        self.trend_p = result.p
        self.trend_significant = result.h

        if result.trend in ("increasing", "decreasing"): ## expected behavior quick exit
            return None 
        
        ### The rest of this is handling the case of result.trend = "no trend" -- untold how common this is 
        # Interesting / problematic case of result.trend = "no trend"
        print("[WARN] Non-monotonic or insignificant trend detected")
        print(f"  Alignment: {self.reference} vs {self.query}")
        print(f"  Trend: {result.trend}")
        print(f"  p-value: {result.p:.4g}")
        print(f"  Significant: {result.h}")
        print(f"  Please see {debug_dir} for diagnostic plots")
        print(f"  Diagnostic reminder - function returns parameter set in this case")
        print(f"  Its probably a good time to get into contact with Trevor =D")

        # Debug plotting -- plot cases where no trend is found. Who knows what's in store for us here. 
        if debug_dir is not None:
            os.makedirs(debug_dir, exist_ok=True)

            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(cumsum)
            ax.set_xlabel("Block index")
            ax.set_ylabel("Cumulative weighted length")
            ax.set_title(
                f"{self.reference} vs {self.query}\n"
                f"Trend: {result.trend}, p={result.p:.3g}"
            )
            ax.grid(True)

            outpath = os.path.join(
                debug_dir,
                f"{self.reference}_vs_{self.query}_trend_debug.png"
            )
            fig.savefig(outpath, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"[INFO] Debug plot written to {outpath}")
        if cumsum[0] <= cumsum[-1]:
            self.trend = "increasing"
        else:
            self.trend = "decreasing"
        self.trend_checked = True
        # ----------------------------------
        # Return diagnostic payload
        # ----------------------------------
        return {
            "reference": self.reference,
            "query": self.query,
            "trend": result.trend,
            "p_value": result.p,
            "significant": result.h
        }
        
    def find_source_sink_nodes(self):
        """
        Determine source and sink blocks for graph traversal.

        Assigns:
            self.source
            self.sink
        """
        if not getattr(self, "trend_checked", False):
            raise RuntimeError(
                "Alignment trend not evaluated; call evaluate_alignment_trend() first"
            )

        if not hasattr(self, "bounding_box"):
            raise RuntimeError(
                "Bounding box not evaluated; call compute_bounding_box() first"
            )
        
        # Choose corners based on trend
        if self.trend == "increasing":
            start = self.bounding_box.bottom_left
            end   = self.bounding_box.top_right
        elif self.trend == "decreasing":
            start = self.bounding_box.top_left
            end   = self.bounding_box.bottom_right
        else:
            raise RuntimeError(
                "Cannot determine source/sink for non-monotonic alignment"
            )

        current_min_start = None
        current_min_end = None
        min_dist_start = float("inf")
        min_dist_end = float("inf")

        for block in self.alignment_blocks:
            d_start = Alignment.x_distance(start, (block.reference_start, block.query_start))
            d_end = Alignment.x_distance(end, (block.reference_end, block.query_end))

            if d_start < min_dist_start:
                min_dist_start = d_start
                current_min_start = block

            if d_end < min_dist_end:
                min_dist_end = d_end
                current_min_end = block

        if current_min_start is None or current_min_end is None:
            raise RuntimeError("Failed to identify source/sink blocks (do you have an alignment with no alignment blocks?)")

        if current_min_start is current_min_end:
            raise RuntimeError("Degenerate traversal: source == sink. Not necessarily fatal but erroring out")

        self.source = current_min_start
        self.sink = current_min_end

class AlignmentBlock:
    def __init__(self, reference, reference_start, reference_end, query_start, query_end, errors, similarity_errors, stop_codons, query , indel_map, index):
        self.reference = reference
        self.query = query
        self.reference_start = int(reference_start)
        self.reference_end = int(reference_end)
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.errors = int(errors)
        self.similarity_errors = int(similarity_errors)
        self.stop_codons = int(stop_codons)
        self.slope = (self.query_end - self.query_start) / (self.reference_end - self.reference_start)
        self.rounded_slope = round(self.slope, 0)
        self.length = self.reference_end - self.reference_start
        self.left_most = min(self.reference_start, self.reference_end)
        self.right_most = max(self.reference_start, self.reference_end)
        self.domain = (self.left_most, self.right_most)
        self.length_as_sole_coverage = 0
        self.y_intercept = self.query_start - (self.slope * self.reference_start) ## b = y - mx
        self.domain_length = self.reference_end - self.reference_start
        self.x_mid = (self.reference_start + self.reference_end) / 2
        self.y_mid = (self.query_start + self.query_end) / 2
        self.z_scaled_y_intercept = None
        self.start = (self.reference_start, self.query_start)
        self.end = (self.reference_end, self.query_end)
        self.left_neighbor = None
        self.right_neighbor = None
        self.shares_domain = []
        self.shares_range = []
        self.part_of_primary_synteny = False
        self.is_repeat = False
        self.repeat_content = 0
        self.index = index
        self.weighted_length = self.slope * self.length
        self.y1_reflection = None
        self.y2_reflection = None
        self.pivot = self.y_mid
        self.reflected = False
        self.indel_map = indel_map
        self.segments = None
        self.segment_tree = None
        
    def __str__(self):
        return f'{self.reference_start} {self.reference_end} {self.query_start} {self.query_end} {self.errors} {self.similarity_errors} {self.stop_codons} {self.slope}'

    def plot_slope(self):
        if self.slope > 0:
            plt.plot([self.reference_start, self.reference_end], [self.query_start, self.query_end], color= 'black')
            #plt.text(self.x_mid, self.y_mid, f'Slope: {self.slope:.4f}', fontsize=5, color='black')
        else:
            plt.plot([self.reference_start, self.reference_end], [self.query_start, self.query_end], color='red')
        
    def plot(self, ax, color='black', text='', title=''):
        if self.is_repeat:
            color = 'green'

        ax.plot(
            [self.reference_start, self.reference_end],
            [self.query_start, self.query_end],
            color=color
        )

        if title:
            ax.set_title(title)

    def plot2(self, color = 'pink', text = '', title = ''):
        plt.plot([self.reference_start, self.reference_end], [self.y1_reflection, self.y2_reflection], color=color)

    def subplot(self, color = 'black', text = '', axis = 0):
        axs[axis].plot([self.reference_start, self.reference_end], [self.query_start, self.query_end], color=color)
        axs[axis].text(self.x_mid, self.y_mid, text, fontsize=10, ha='right')

    def build_segments(self):
        """
        Derive ungapped syntenic segments from delta indel map.
        Effectively deprecated, kept in case I explore constituent line segments in the future
        """
        if self.segments is not None:
            return

        ref_pos = self.reference_start
        qry_pos = self.query_start
        self.segments = []

        for d in self.indel_map:
            if d == 0:
                break

            seg_len = abs(d)
            if seg_len > 1:
                self.segments.append(
                    LineSegment(
                        ref_start=ref_pos,
                        ref_end=ref_pos + seg_len,
                        query_start=qry_pos,
                        query_end=qry_pos + seg_len
                    )
                )

            ref_pos += seg_len
            qry_pos += seg_len

            if d > 0:
                ref_pos += 1
            else:
                qry_pos += 1

    def build_interval_tree(self):
        """Build reference-space interval tree over segments."""
        if self.segment_tree is not None:
            return

        if self.segments is None:
            raise RuntimeError("build_segments() must be called first")

        self.segment_tree = IntervalTree()
        for seg in self.segments:
            self.segment_tree[seg.ref_start:seg.ref_end] = seg

@dataclass
class LineSegment:
    '''Effectively deprecated, kept in case I explore each MUM's constituent line segments in the future'''
    ref_start: int
    query_start: int
    ref_end: int
    query_end: int

    def __str__(self):
        return f"{self.ref_start} {self.query_start} {self.ref_end} {self.query_end}"
    
    def plot(self, ax = None, color = None, linewidth = 1.5):
        from matplotlib import pyplot as plt
        # Create axis if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))

        # Default color if none specified
        if color is None:
            color = "blue"

        # Extract coordinates
        x = [self.ref_start, self.ref_end]
        y = [self.query_start, self.query_end]

        # Plot the segment
        ax.plot(x, y, color=color, linewidth=linewidth)

        return ax

class DeltaParser:
    @staticmethod
    def parse_file(file_path, threads = 1):
        alignments = []
        try:
            it = peekable(open(file_path)) # allows us to "peek" next index w/o consuming
            next(it) # consume first two lines (file names and software info)
            next(it)
            aln = None
            for line in it:
                line = line.strip()
                if line.startswith('>'):
                    index = 0
                    if aln:
                        alignments.append(aln)
                    parts = line.strip('>').split()
                    aln = Alignment(parts[0], parts[1], parts[2], parts[3])

                elif len(line.split()) == 7:  # MUM header
                    indel_map = []
                    while True:
                        nxt = next(it).strip()
                        val = int(nxt)
                        indel_map.append(val)
                        if val == 0:
                            break

                    aln.add_alignment_block(line.split(), indel_map, index)
                    index += 1

                else:
                    print("There has been an error, potentially a malformed delta file")
                    sys.exit(1)
    
            if aln:  # Append the last alignment if exists
                alignments.append(aln)

            return alignments
        
        except FileNotFoundError:
            print(f'File {file_path} not found. Aborting.')
            sys.exit(1)

def select_primary_alignments(alignments):
    ''' We can move this to a generic core module if needed later
    '''
    by_ref = defaultdict(list)
    for aln in alignments:
        by_ref[aln.reference].append(aln)
    primary = {}

    for ref, alns in by_ref.items():
        for aln in alns:
            aln.calculate_coverage(aln.reference_length)

        primary[ref] = max(alns, key=lambda a: a.covered_sites)

    return primary

def plot_alignment_blocks(aln, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))

    for block in aln.alignment_blocks:
        block.plot(ax=ax)

    ax.set_xlabel(f"{aln.reference} (len={aln.reference_length})")
    ax.set_ylabel(f"{aln.query} (len={aln.query_length})")

    return ax

def plot_bounding_box(aln, ax, color="red", linestyle="--", linewidth=1.5):
    """
    Overlay the alignment bounding box on an existing dotplot axis.
    """
    if not hasattr(aln, "bounding_box"):
        raise RuntimeError("Bounding box not computed; call compute_bounding_box() first")

    bb = aln.bounding_box

    xs = [
        bb.bottom_left[0],
        bb.top_left[0],
        bb.top_right[0],
        bb.bottom_right[0],
        bb.bottom_left[0],  # close the box
    ]
    ys = [
        bb.bottom_left[1],
        bb.top_left[1],
        bb.top_right[1],
        bb.bottom_right[1],
        bb.bottom_left[1],
    ]

    ax.plot(xs, ys, color=color, linestyle=linestyle, linewidth=linewidth)


def extract_alignment_block_fastas(aln, reffasta, qryfasta, outdir=None):
    """
    Extract reference/query FASTAs for each alignment block.

    Returns a list of dicts with block metadata.
    """

    records = []

    for index, block in enumerate(aln.alignment_blocks):
        ref_chr = reffasta[block.reference]
        qry_chr = qryfasta[block.query]

        ref_start = min(block.reference_start, block.reference_end)
        ref_end   = max(block.reference_start, block.reference_end)
        ref_seq = ref_chr[ref_start:ref_end]

        qry_start = min(block.query_start, block.query_end)
        qry_end   = max(block.query_start, block.query_end)
        qry_raw = qry_chr[qry_start:qry_end]

        qry_seq = (
            qry_raw.reverse.complement
            if block.slope < 0
            else qry_raw
        )

        if outdir is not None:
            outname = f"{outdir}/{block.reference}_{index}.fasta"
            with open(outname, "w") as out:
                out.write(">reference\n")
                out.write(str(ref_seq) + "\n")
                out.write(">query\n")
                out.write(str(qry_seq) + "\n")

        records.append({
            "block_index": index,
            "reference": block.reference,
            "ref_start": ref_start,
            "ref_end": ref_end,
            "query": block.query,
            "qry_start": qry_start,
            "qry_end": qry_end,
            "length": ref_end - ref_start,
            "slope": block.slope,
        })

    return records

def summarize_block_sizes(block_records):
    import numpy as np

    lengths = np.array([r["length"] for r in block_records])

    return {
        "n_blocks": len(lengths),
        "median": np.median(lengths),
        "mean": np.mean(lengths),
        "p25": np.percentile(lengths, 25),
        "p75": np.percentile(lengths, 75),
        "max": np.max(lengths),
        "lengths": lengths,
    }

def plot_block_size_distribution(lengths):
    import numpy as np
    import matplotlib.pyplot as plt

    lengths = np.sort(lengths)
    y = np.arange(1, len(lengths) + 1) / len(lengths)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    # Histogram
    ax[0].hist(lengths, bins=50)
    ax[0].set_xlabel("Block size (bp)")
    ax[0].set_ylabel("Count")

    # ECDF
    ax[1].plot(lengths, y)
    ax[1].set_xlabel("Block size (bp)")
    ax[1].set_ylabel("Fraction ≤ size")
    ax[1].set_xscale("log")

    return ax


def main():
    # ----------------------------
    # CLI facing main 
    # ----------------------------
    parser = argparse.ArgumentParser(
        description="Parse a MUMmer delta file and plot primary synteny alignments"
    )
    parser.add_argument(
        "delta",
        help="Input .delta file"
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Output directory for plots (default: current directory)"
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display plots interactively (default: save only)"
    )

    args = parser.parse_args()

    delta = args.delta
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    # ----------------------------
    # Parse delta
    # ----------------------------
    alignments = DeltaParser.parse_file(delta)
    primary = select_primary_alignments(alignments)

    if not primary:
        print("[WARN] No primary alignments found")
        sys.exit(0)

    # ----------------------------
    # Plot each primary alignment
    # ----------------------------
    for key, aln in primary.items():
        print(f"[INFO] Processing alignment: {key}")

        aln.build_synteny()

        fig, ax = plt.subplots(figsize=(10, 6))

        # Base alignment blocks
        plot_alignment_blocks(aln, ax=ax)
        plot_bounding_box(aln, ax=ax)

        # Overlay primary synteny blocks (sanity / emphasis)
        for block in aln.primary_synteny_blocks:
            block.plot(ax=ax, color="blue")

        ax.set_xlabel(f"{aln.reference} ({aln.reference_length:,} bp)")
        ax.set_ylabel(f"{aln.query} ({aln.query_length:,} bp)")
        ax.set_title(f"{aln.reference} vs {aln.query}")

        outpath = os.path.join(outdir, f"{key}_synteny.png")
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
        print(f"[✔] Wrote {outpath}")

        if args.show:
            plt.show()

        plt.close(fig)



if __name__ == '__main__':
    main()
