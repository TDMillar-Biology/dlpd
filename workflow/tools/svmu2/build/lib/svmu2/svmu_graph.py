'''
TDMILLAR SVMU2.0
Goal: Detect Structural Variation in Genomes via Whole Genome Alignment
Data: Whole Genome Alignments produced with MUMMER (delta format)
Strategy: Create a densly connected directed graph, Use Djikstras algorithm to find best traversal
Output: VCF File indicating the coordinates of SVs, could include BED format and others
Contact: Tmillar32@gmail.com

Jan 15 2025 revisions
'''

from delta_parser import DeltaParser
#from .dotplot_line_segment import DotPlotLineSegment
#from .parse_repeat_masker import ParseRepeatMasker
import matplotlib.pyplot as plt
from math import sqrt
import plotly.io as pio
import plotly.tools as pt
import math
import plotly.graph_objects as go
from intervaltree import IntervalTree
import argparse
import os
from statistics import mode
import pysam
from datetime import datetime
import sys
import numpy as np
from collections import defaultdict
import pymannkendall as mk
import pdb
import networkx as nx
from math import sqrt
from sortedcontainers import SortedList
from statistics import mode

def plot_static_dotplot(best_delta, new_SVs):
    log(f"Plotting static dotplot for {best_delta.reference}")
    fig, ax = plt.subplots()

    for sv in new_SVs:
        color = 'yellow' if sv.sv_type == 'BND' else 'orange'
        sv.plot(color=color)

    for block in best_delta.alignment_blocks:
        block.plot(color='blue')

    for block in best_delta.alignment_blocks:
        if not block.is_repeat and not block.part_of_primary_synteny:
            block.plot(color='pink')

    for block in best_delta.alignment_blocks:
        if block.part_of_primary_synteny:
            block.plot(color='blue')

    fig.savefig(f'{best_delta.reference}_dotplot.png')
    log(f"Saved plot: {best_delta.reference}_dotplot.png")
    plt.close(fig)

def plot_interactive_dotplot(best_delta, new_SVs=None, x_marker=None):
    log(f"Plotting interactive dotplot for {best_delta.reference}")
    plotly_data = []

    # Add SVs as scatter lines
    if new_SVs:
        for sv in new_SVs:
            color = 'yellow' if sv.sv_type == 'BND' else 'orange'
            trace = go.Scatter(
                x=[sv.reference_start, sv.reference_end],
                y=[sv.query_start, sv.query_end],
                mode='lines',
                line=dict(color=color),
                hovertext=f"SV type: {sv.sv_type}",
                hoverinfo='text'
            )
            plotly_data.append(trace)

    for block in best_delta.alignment_blocks:
        # Always plot the original block
        block_color = (
            'grey' if (not block.is_repeat and not block.part_of_primary_synteny) else
            'blue' if block.part_of_primary_synteny else
            'black'
        )
        hover_text = (
            f"Index: {block.index}<br>"
            f"Ref: {block.reference_start}-{block.reference_end}<br>"
            f"Query: {block.query_start}-{block.query_end}<br>"
            f"Repeat: {block.is_repeat}<br>"
            f"Primary Synteny: {block.part_of_primary_synteny}<br>"
            f"Reflected: {getattr(block, 'reflected', False)}"
        )
        orig_trace = go.Scatter(
            x=[block.reference_start, block.reference_end],
            y=[block.query_start, block.query_end],
            mode='lines',
            line=dict(color=block_color),
            hovertext=hover_text,
            hoverinfo='text',
        )
        plotly_data.append(orig_trace)

        # If reflected, also plot the reflected coordinates
        if block.reflected:
            reflect_trace = go.Scatter(
                x=[block.reference_start, block.reference_end],
                y=[block.y1_reflection, block.y2_reflection],
                mode='lines',
                line=dict(color='green'),
                hovertext = (
                    f"Index: {block.index}<br>"
                    f"Ref: {block.reference_start}-{block.reference_end}<br>"
                    f"Query: {block.y1_reflection}-{block.y2_reflection}<br>"
                    f"Repeat: {block.is_repeat}<br>"
                    f"Primary Synteny: {block.part_of_primary_synteny}<br>"
                    f"Reflected: {getattr(block, 'reflected', False)}"
                ),
                hoverinfo='text',
            )
            print('working')
            plotly_data.append(reflect_trace)

    # Add vertical line at x_marker (e.g., x_target for regression line)
    if x_marker is not None:
        miny = min(b.query_start for b in best_delta.alignment_blocks)
        maxy = max(b.query_end for b in best_delta.alignment_blocks)
        line_trace = go.Scatter(
            x=[x_marker, x_marker],
            y=[miny, maxy],
            mode='lines',
            line=dict(color='red', dash='dash'),
            name='X Marker',
            hoverinfo='skip'
        )
        plotly_data.append(line_trace)

    layout = go.Layout(
        title=f"Interactive Dotplot: {best_delta.reference}",
        xaxis=dict(title='Reference'),
        yaxis=dict(title='Query'),
        showlegend=False
    )
    
    plotly_fig = go.Figure(data=plotly_data, layout=layout)
    pio.write_html(plotly_fig, file=f'{best_delta.reference}_interactive.html', auto_open=True)
    log(f"Saved interactive plot: {best_delta.reference}_interactive.html")

def plot_bw_inset(best_delta):
    log(f"Generating black and white plot for {best_delta.reference}")
    fig, ax = plt.subplots()

    for block in best_delta.alignment_blocks:
        block.plot(color='black')

    ax.set_xlabel(f'{best_delta.query} Reference')
    ax.set_ylabel(f'{best_delta.reference} Query')
    ax.set_title(r"$\it{D.\ melanogaster}$ " + f"{best_delta.reference}")

    fig.savefig(f'{best_delta.reference}_dotplot_bw_inset.png')
    log(f"Saved black and white plot: {best_delta.reference}_dotplot_bw_inset.png")
    plt.close(fig)

class AlignmentSubset:
    def __init__(self, blocks):
        self.alignment_blocks = blocks
        self.query_length = max(b.query_end for b in blocks)
        self.reference_length = max(b.reference_end for b in blocks)
        self.slope = None

def calculate_distance(aln1, aln2):
    ''' Calulate the distance between the end of aln1 and the beginning of aln2
    '''
    x1, y1 = aln1.end
    x2, y2 = aln2.start
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def euclidean_ordered_pairs(a, b):
    ''' a and b are x,y tupples
    '''
    return math.sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2)

def pairwise_euclidean_distance(object, collection):
    distances = {}
    for obj in collection:
        distance = calculate_distance(object, obj)
        distances[distance] = obj
    return distances

def find_theta(obj1, obj2):
    ''' 
    Calculate the angle theta (in degrees) in a polar coordinate plane 
    that defines the direction of the line segment from obj1 to obj2.
    '''
    # Calculate change in x and y
    delta_x = obj2.reference_start - obj1.reference_end
    delta_y = obj2.query_start - obj1.query_end

    # Use atan2 for proper angle computation across all quadrants
    angle_rad = math.atan2(delta_y, delta_x) # div by zero protection returns (-pi,pi) 
    angle_deg = math.degrees(angle_rad)
    
    return angle_deg

def is_small(x, threshold = 20): #20 nucleotides is what mummer is capable of detecting per publication
    return abs(x) < threshold

def classify_sv(segment, domain_tree, range_tree, write_bnds=False):
    """
    Classifies the SV type using vector geometry and domain/range interval trees.
    Returns a list of DotPlotLineSegments (primary + optional breakends) sharing an event ID.
    """
    x_component = segment.reference_end - segment.reference_start
    y_component = segment.query_end - segment.query_start


    domain_hits = domain_search(segment, domain_tree)
    segment.domain_partners = domain_hits

    range_hits = range_search(segment, range_tree)
    segment.range_partners = range_hits

    # Assign a unique event ID for this group of calls
    event_id = f'event_{segment.chrom}_{segment.index:04}'
    segment.event_ID = event_id

    sv_records = [segment]  # list of DotPlotLineSegment to return

    # Case 1: Too small — skip writing
    if is_small(x_component) and is_small(y_component):
        segment.sv_type = 'UNDETERMINED'
        return sv_records

    # Case 3: Large X, small Y — DEL or DUP (reference-space event)
    if is_small(y_component):
        if abs(segment.theta) < 45:
            segment.sv_type = 'DEL'
            # Search range for potential relocated sequence
            if write_bnds:
                for start, stop, hit in range_hits:
                    bnd = DotPlotLineSegment(
                        chrom=segment.chrom,
                        reference_start=hit.reference_start,
                        reference_end=hit.reference_end,
                        query_start=hit.query_start,
                        query_end=hit.query_end,
                        sv_type='BND',
                        theta=segment.theta,
                        index=segment.index
                    )
                    bnd.event_ID = event_id
                    sv_records.append(bnd)
        elif abs(abs(segment.theta) - 180) < 45:
            if len(domain_hits) > 2:  # Check for tandem duplication
                segment.sv_type = 'DUP_TANDEM'
            else:
                segment.sv_type = 'DUP'
        else:
            segment.sv_type = 'COMPLEX'
        return sv_records

    # Case 4: Small X, large Y — INS or tandem DUP (query-space event)
    if is_small(x_component):
        if abs(segment.theta - 90) < 45:
            segment.sv_type = 'INS'
            # Search domain for potential insertion source
            if write_bnds:
                for start, stop, hit in domain_hits:
                    bnd = DotPlotLineSegment(
                        chrom=segment.chrom,
                        reference_start=hit.reference_start,
                        reference_end=hit.reference_end,
                        query_start=hit.query_start,
                        query_end=hit.query_end,
                        sv_type='BND',
                        theta=segment.theta,
                        index=segment.index
                    )
                    bnd.event_ID = event_id
                    sv_records.append(bnd)
        elif abs(segment.theta + 90) < 45:
            if len(domain_hits) > 2:  # Check for tandem duplication
                segment.sv_type = 'DUP_TANDEM'
            else:
                segment.sv_type = 'DUP'
        else:
            segment.sv_type = 'COMPLEX'
        return sv_records

    # Case 5: Not aligned with x or y axes → complex this currently works as expected
    if abs(segment.theta) > 90 and abs(segment.theta) < 180:
        if len(domain_hits) > 2:  # Check for tandem duplication
            segment.sv_type = 'DUP_TANDEM'
        else:
            segment.sv_type = 'DUP'
    else:
        segment.sv_type = 'COMPLEX'
        
    if write_bnds:
        for start, stop, hit in domain_hits:
            bnd = DotPlotLineSegment(
                chrom=segment.chrom,
                reference_start=hit.reference_start,
                reference_end=hit.reference_end,
                query_start=hit.query_start,
                query_end=hit.query_end,
                sv_type='BND',
                theta=segment.theta,
                index=segment.index
            )
            bnd.event_ID = event_id
            sv_records.append(bnd)

        for start, stop, hit in range_hits:
            bnd = DotPlotLineSegment(
                chrom=segment.chrom,
                reference_start=hit.reference_start,
                reference_end=hit.reference_end,
                query_start=hit.query_start,
                query_end=hit.query_end,
                sv_type='BND',
                theta=segment.theta,
                index=segment.index
            )
            bnd.event_ID = event_id
            sv_records.append(bnd)
    return sv_records

def write_vcf(SVs, output_path="output.vcf", sample='SAMPLE'):
    import pysam
    header = pysam.VariantHeader()
    header.add_meta('INFO', items=[('ID', 'SVTYPE'), ('Number', '1'), ('Type', 'String'), ('Description', 'Type of structural variant')])
    header.add_meta('INFO', items=[('ID', 'END'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'End position of the SV')])
    header.add_meta('INFO', items=[('ID', 'SVLEN'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Length of the SV')])
    header.add_meta('INFO', items=[('ID', 'THETA'), ('Number', '1'), ('Type', 'Float'), ('Description', 'Angle between segments in degrees')])
    header.add_meta('INFO', items=[('ID', 'IMPRECISE'), ('Number', '0'), ('Type', 'Flag'), ('Description', 'Imprecise structural variation')])
    header.add_meta('INFO', items=[('ID', 'DUP_TANDEM'), ('Number', '0'), ('Type', 'Flag'), ('Description', 'Tandem duplication')])
    header.add_meta('INFO', items=[('ID', 'CN'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Estimated copy number in query')])
    header.add_meta('INFO', items=[('ID', 'CNL'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Estimated copy number in reference')])
    header.add_meta('INFO', items=[('ID', 'EVENT'),('Number', '1'),('Type', 'String'),('Description', 'ID of the complex event this variant is part of')])
    
    header.add_sample(sample)

    contigs = sorted(set(SV.chrom for SV in SVs)) # get each unique contig
    for contig in contigs:
        header.contigs.add(contig) ## add them to the vcf header

    vcf_out = pysam.VariantFile(output_path, 'w', header=header)

    for i, SV in enumerate(SVs): #########################
        if SV.theta == 45 and not (SV.range_partners or SV.domain_partners):
            continue
        chrom = SV.chrom

        if SV.theta < -90 or SV.theta > 90:
            pos = SV.reference_end
            end = SV.reference_start
        else:
            pos = SV.reference_start
            end = SV.reference_end

        length = end - pos if SV.sv_type != "INS" else SV.query_end - SV.query_start
        ref_base = 'N'
        alt = f"<{SV.sv_type if SV.sv_type != 'DUP_REF' else 'DEL'}>"

        info_dict = {
            'SVTYPE': SV.sv_type if SV.sv_type != 'DUP_REF' else 'DEL',
            'END': end,
            'SVLEN': length,
            'EVENT': SV.event_ID,
            'THETA': SV.theta,
            'IMPRECISE': True
        }

        # Annotate tandem duplication flags
        if SV.sv_type in ('DUP', 'DUP_REF'):
            info_dict['DUP_TANDEM'] = True
            info_dict['CN'] = 2 if SV.sv_type == 'DUP' else 1
            info_dict['CNL'] = 1 if SV.sv_type == 'DUP' else 2

        record = vcf_out.new_record(
            contig=chrom,
            start=pos,
            stop=end,
            alleles=(ref_base, alt),
            id=f'sv_{i}',
            info=info_dict
        )
        vcf_out.write(record)

    vcf_out.close()

def dijkstra_traversal(blocks, source, sink, max_jump=100_000, weights='euclidean', return_graph=False):
    """
    Perform a Dijkstra traversal on alignment blocks with optional orientation-aware weights.

    Parameters:
        blocks (list): List of alignment block objects.
        source (int): Index of the source node in blocks.
        sink (int): Index of the sink node in blocks.
        max_jump (int): Max reference distance to consider for edges.
        weights (str): One of 'euclidean', 'X', or 'Y'.
        return_graph (bool): Whether to include the full graph in the output.

    Returns:
        dict: Contains source, sink, path (indices), total_distance, path_blocks, and optionally the graph.
    """
    def smart_euclidean_distance(a, b):

        return orientation_aware_euclidean_distance(a, b)

    def euclidean_distance(x1, y1, x2, y2):
        return sqrt((x2 - x1)**2 + (y2 - y1)**2)

    def orientation_aware_euclidean_distance(a, b):
        scores = []
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_start))
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_end))
        scores.append(euclidean_distance(a.reference_end, a.query_start, b.reference_start, b.query_start))
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_end, b.query_start))
        if getattr(a, "reflected", False):
            scores.append(euclidean_distance(a.reference_end, a.y2_reflection, b.reference_start, b.query_start)) # reflection
            #scores.append(euclidean_distance(a.reference_start, a.y1_reflection, b.reference_start, b.query_start)) # reflection
        if getattr(b, "reflected", False):
            scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.y1_reflection)) # reflection
            #scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_end, b.y2_reflection)) # reflection
        min_index = min(enumerate(scores), key=lambda x: x[1])[0]
        with open('distances.txt', 'a') as out:
            out.write(f'{a.index}, {b.index}, {min(scores)},{min_index}\n')

        return min(scores)
    
    def xdistance(a, b):
        return abs(b.reference_start - a.reference_end)

    def ydistance(a, b):
        return abs(b.query_start - a.query_end)

    def get_reference_neighbors(ref_index, end_coord, window):
        left = ref_index.bisect_left((end_coord - window, -1))
        right = ref_index.bisect_right((end_coord + window, float('inf')))
        return ref_index[left:right]

    # Build reference-based index for neighbor search
    reference_index = SortedList([(block.reference_start, i, block) for i, block in enumerate(blocks)])

    G = nx.DiGraph()
    for i, block in enumerate(blocks):
        G.add_node(i, block=block)

    weight_fn = {
        'euclidean': smart_euclidean_distance,
        'X': xdistance,
        'Y': ydistance
    }.get(weights)

    if weight_fn is None:
        raise ValueError(f"Invalid weight designation for Dijkstra traversal: {weights}")

    # Expand search window gradually if needed
    for factor in [1, 3, 5]:
        current_jump = max_jump * factor

        for i, a in enumerate(blocks):
            neighbors = get_reference_neighbors(reference_index, a.reference_end, current_jump)
            for _, j, b in neighbors:
                if i == j:
                    continue
                distance = weight_fn(a, b)
                G.add_edge(i, j, weight=distance)

        try:
            path = nx.dijkstra_path(G, source, sink, weight='weight')
            total_distance = sum(G[path[i]][path[i + 1]]['weight'] for i in range(len(path) - 1))
        except nx.NetworkXNoPath:
            path = []
            total_distance = float('inf')

        if path:
            break

    result = {
        'source': source,
        'sink': sink,
        'path': path,
        'total_distance': total_distance,
        'path_blocks': [G.nodes[i]['block'] for i in path]
    }

    if return_graph:
        result['graph'] = G

    return result

def extract_dotplot_segments_from_path(block_path, slope, reference_name, domain_tree, range_tree, write_bnds=False):
    all_segments = []
    for i in range(len(block_path) - 1):
        a, b = block_path[i], block_path[i + 1]
        theta = find_theta(a, b)
        if a.rounded_slope == b.rounded_slope == slope:
            segment = DotPlotLineSegment(
                chrom=reference_name,
                reference_start=a.reference_end,
                reference_end=b.reference_start,
                query_start=a.query_end,
                query_end=b.query_start,
                sv_type=None,
                theta=theta,
                index=i
            )
            if write_bnds:
                segments = classify_sv(segment, domain_tree, range_tree, write_bnds=True)
            else:
                segments = classify_sv(segment, domain_tree, range_tree)

            all_segments.extend(segments)
        elif a.rounded_slope == slope and a.rounded_slope != b.rounded_slope:
            # enter inversion from left DO NOTHING 
            pass

    return all_segments
    
def create_domain_range_trees(alignment_blocks, min_overlap_buffer=25):
    '''
    Accept Aligment blocks as iterable
    ensures start < end as interval tree requires
    returns tuple of objects being (domain_tree, range_tree)
    This implementation of an interval tree allows me to search for overlapping segments in log n time
    '''
    domain_tree = IntervalTree()
    range_tree = IntervalTree()

    for block in alignment_blocks:
        # Reference (domain) interval
        ref_start = min(block.reference_start, block.reference_end)
        ref_end = max(block.reference_start, block.reference_end)

        if (ref_end - ref_start) > 2 * min_overlap_buffer:
            ref_start += min_overlap_buffer
            ref_end   -= min_overlap_buffer
            domain_tree[ref_start:ref_end] = block
        else:
            #print(f"Skipping short reference block: {block}")
            pass

        # Query (range) interval
        qry_start = min(block.query_start, block.query_end)
        qry_end = max(block.query_start, block.query_end)

        if (qry_end - qry_start) > 2 * min_overlap_buffer:
            qry_start += min_overlap_buffer
            qry_end   -= min_overlap_buffer
            range_tree[qry_start:qry_end] = block
        else:
            #print(f"Skipping short query block: {block}")
            pass

    return domain_tree, range_tree

def create_repeat_tree(repeats, chrom, buffer = 5):
    ''' 
    Accepts output repeats from ParseRepeatMasker as input
    returns interval tree of (start, stop, repeatobject)    
    '''
    repeat_tree = IntervalTree()
    filtered_repeats = [r for r in repeats if r.chrom == chrom]

    for repeat in filtered_repeats:
        qry_start = int(repeat.start) - buffer
        qry_end = int(repeat.end) + buffer

        # Validate
        if qry_end > qry_start:
            repeat_tree[qry_start:qry_end] = repeat
        else:
            print(f"Skipping invalid or shrunken repeat: {repeat.start}-{repeat.end} → {qry_start}-{qry_end}")

    return repeat_tree

def domain_search(dotplot_segment, tree, fraction_threshold=0.50):
    """
    Returns a list of (overlap_length, block) for reference-overlapping blocks,
    where the overlap covers at least fraction_threshold of dotplot_segment's reference span.
    """
    seg_start = min(dotplot_segment.reference_start, dotplot_segment.reference_end)
    seg_end   = max(dotplot_segment.reference_start, dotplot_segment.reference_end)
    seg_len = seg_end - seg_start
    if seg_len <= 0:
        return []

    results = []
    for interval in tree[seg_start:seg_end]:
        overlap_start = max(seg_start, interval.begin)
        overlap_end = min(seg_end, interval.end)
        overlap_length = max(0, overlap_end - overlap_start)
        if overlap_length / seg_len >= fraction_threshold:
            results.append((overlap_length, interval.data))  # .data is a block
    return results

def range_search(dotplot_segment, tree, fraction_threshold=0.50):
    """
    Returns a list of (overlap_length, block) for query-overlapping blocks,
    where the overlap covers at least fraction_threshold of dotplot_segment's query span.
    """
    seg_start = min(dotplot_segment.query_start, dotplot_segment.query_end)
    seg_end   = max(dotplot_segment.query_start, dotplot_segment.query_end)
    seg_len = seg_end - seg_start
    if seg_len <= 0:
        return []

    results = []
    for interval in tree[seg_start:seg_end]:
        overlap_start = max(seg_start, interval.begin)
        overlap_end = min(seg_end, interval.end)
        overlap_length = max(0, overlap_end - overlap_start)
        if overlap_length / seg_len >= fraction_threshold:
            results.append((overlap_length, interval.data))  # .data is a block
    return results

def repeat_search(alignment_block, tree):
    return tree[alignment_block.query_start:alignment_block.query_end]

def annotate_repeats_with_fraction(blocks, repeat_tree, fraction_threshold=0.90):
    for block in blocks:
        query_start = min(block.query_start, block.query_end)
        query_end = max(block.query_start, block.query_end)
        block_length = query_end - query_start
        if block_length <= 0:
            continue  # Skip malformed blocks

        # Get overlapping repeats
        overlaps = repeat_tree[query_start:query_end]
        if not overlaps:
            continue  # No repeats

        # Build a new temporary interval tree to merge overlaps
        overlap_tree = IntervalTree()
        for interval in overlaps:
            repeat_start = max(query_start, interval.begin)
            repeat_end = min(query_end, interval.end)
            if repeat_end > repeat_start:
                overlap_tree[repeat_start:repeat_end] = True  # placeholder

        # Compute total overlap length from merged intervals
        overlap_tree.merge_overlaps()
        overlap_bases = sum(iv.end - iv.begin for iv in overlap_tree)

        if (overlap_bases / block_length) >= fraction_threshold:
            block.is_repeat = True
            block.repeat_content = (overlap_bases / block_length)
            print(f"Marked as repeat: {query_start}-{query_end} | Overlap: {overlap_bases}/{block_length:.1f}")  
  
def log(msg: str, level: str = "INFO"):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [svmu] [{level}] - {msg}")

def score_path(path_blocks, reflen):
    if reflen == 0:
        return 0
    total_length = sum(block.length for block in path_blocks)
    return total_length / reflen

def recursive_trend_segment(aln, depth=0):
    trace = [b.weighted_length for b in aln.alignment_blocks]
    cumsum = np.cumsum(trace)
    result = mk.hamed_rao_modification_test(cumsum)
    print(depth)
    if result.h:
        if result.trend == 'increasing':
            aln.slope = 1
        elif result.trend == 'decreasing':
            aln.slope = -1
        return[aln]

    # Otherwise, break into parts and recurse
    min_index = np.argmin(cumsum)
    max_index = np.argmax(cumsum)
    start, end = sorted([min_index, max_index])
    print(start, end)
    segments = []
    if start > 1:
        sub_aln1 = AlignmentSubset(aln.alignment_blocks[:start])
        segments.extend(recursive_trend_segment(sub_aln1, depth+1))
    sub_aln2 = AlignmentSubset(aln.alignment_blocks[start:end+1])
    segments.extend(recursive_trend_segment(sub_aln2, depth+1))
    if end+1 < len(aln.alignment_blocks):
        sub_aln3 = AlignmentSubset(aln.alignment_blocks[end+1:])
        segments.extend(recursive_trend_segment(sub_aln3, depth+1))

    return segments

def reflect_alignment_blocks(alignments, trend):
    """
    Reflects alignment blocks across their assigned pivot (x-axis reflection).
    Assumes that each block needing reflection has a `.pivot` attribute set.
    Updates the reference coordinates in-place.
    """
    for aln in alignments:
        for block in aln.alignment_blocks:
            if aln.slope != trend:
                point = find_xeqc_point_intersection(block)
                print(point)
                reflect_line_segment(block, point)

def find_xeqc_point_intersection(block):
    y2 = ((block.pivot - block.reference_start) * block.slope) + block.query_start
    return (block.pivot, y2)

def reflect_line_segment(block, point):
    slope = -1/block.slope
    x1, y1 = point
    block.y2_reflection = ((block.reference_end - x1) * slope) + y1
    block.y1_reflection = ((block.reference_start - x1) * slope) + y1
    block.reflected = True

def find_pivots(alignments, trend):
    for aln in alignments:
        if aln.slope != trend:  # This is an inversion
            minx = min(block.reference_start for block in aln.alignment_blocks)
            maxx = max(block.reference_start for block in aln.alignment_blocks)
            midx = (maxx + minx) / 2
            print(minx, maxx, midx)
            for block in aln.alignment_blocks:
                block.pivot = midx

def compute_primary_synteny(alignment, max_jump=100_000, weights="euclidean"):
    blocks = alignment.alignment_blocks

    try:
        source_idx = blocks.index(alignment.source)
        sink_idx   = blocks.index(alignment.sink)
    except ValueError:
        raise RuntimeError("Source or sink block not found in alignment blocks list")

    traversal = dijkstra_traversal(
        blocks,
        source_idx,
        sink_idx,
        max_jump=max_jump,
        weights=weights
    )

    alignment.set_primary_synteny_blocks(traversal["path_blocks"])
    alignment.build_primary_synteny_tree()


def main():
    parser = argparse.ArgumentParser(prog="svmu2", description="Structural Variants from MUmer")
    parser.add_argument('--delta', '-d', required=True, help='Path to .delta alignment file')
    parser.add_argument('--plot', '-p', action='store_true', help='Generate static dotplot(s)')
    parser.add_argument('--call', '-c', action='store_true', help='Call SVs and write to VCF')
    parser.add_argument('--interactive', '-i', action='store_true', help='Launch interactive plot viewer')
    parser.add_argument('--scaffold', '-sc', action='store_true', help='Perform reference-based scaffolding')
    parser.add_argument('--output', '-o', help='Output VCF file (required if --call is used)')
    parser.add_argument('--reference', '-r', help='Reference FASTA (required if --scaffold)')
    parser.add_argument('--query', '-q', help='Query FASTA (required if --scaffold)')
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    parser.add_argument('--sample', '-s', default='SAMPLE', help='Sample name for VCF output')
    parser.add_argument('--repeatmasker', '-rm', help='Path to RepeatMasker .out file for repeat annotation')
    parser.add_argument('--max_jump', '-mj', type=int, default=100_000, help='Maximum distance to consider for connecting blocks (see README before adjusting)')
    parser.add_argument('--write-bnds', '-wb', action='store_true', help='Also output BND breakends for each SV (default: False, niche use case, not recommended for most users)')
    parser.add_argument('--plotbw', '-bw', action='store_true', help='Generate static dotplot(s) in black and white')
    parser.add_argument('--reference_contigs', '-rc', action='store_true', help='Print reference contigs found in delta file')
    parser.add_argument("--canonical_contigs", help="File or comma-separated list of canonical contig names.")
    args = parser.parse_args()

    if not (args.plot or args.call or args.interactive or args.scaffold or args.plotbw):
        log("No action specified. Use at least one of --plot, --call, --interactive, or --scaffold. or --plotbw", level="ERROR")
        sys.exit(1)
    if args.call and not args.output:
        log("--call requires --output", level="ERROR")
        sys.exit(1)
    if args.scaffold and (not args.reference or not args.query):
        log("--scaffold requires --reference and --query", level="ERROR")
        sys.exit(1)
    if not os.path.exists(args.delta):
        log(f"Delta file not found: {args.delta}", level="ERROR")
        sys.exit(1)
    if args.repeatmasker and not os.path.exists(args.repeatmasker):
        log(f"RepeatMasker file not found: {args.repeatmasker}", level="ERROR")
        sys.exit(1)
    if args.repeatmasker:
        repeats = ParseRepeatMasker.parse_repeat_masker_out(args.repeatmasker)
    else:
        repeats = None
    # Read contig names from file or comma-separated string
    if args.canonical_contigs:
        if "," in args.canonical_contigs:
            canonical_contigs = [x.strip() for x in args.canonical_contigs.split(",")]
        else:
            with open(args.canonical_contigs) as f:
                canonical_contigs = [line.strip() for line in f if line.strip()]
    else:
        log("No canonical contigs provided. Using default for drosophila: [2L, 2R, 3L, 3R, 4, X, Y]", level="WARNING")
        canonical_contigs = ['2L', '2R', '3L', '3R', '4', 'X', 'Y'] # Could set to default for humans: [str(x) for x in range(1,23)] + ['X','Y']

    log("Parsing delta file...")
    deltas = DeltaParser.parse_file(args.delta, threads=args.threads)
    #reference_contigs = sorted(set(delta.reference for delta in deltas))
    #log(f"Delta contains the following reference contigs: {', '.join(reference_contigs)}")

    ref_to_deltas = defaultdict(list)
    for delta in deltas:
        ref_to_deltas[delta.reference].append(delta) # Find all deltas which have a common reference

    for ref in ref_to_deltas:
        ref_to_deltas[ref].sort(key=lambda d: d.covered_sites, reverse=True) # sort those deltas by covered sites 

    SVs = []
    for i, contig in enumerate(canonical_contigs[:7]):
        log(f"Processing alignment {contig} ({i+1}/{len(canonical_contigs)})")
        if ref_to_deltas[contig]:
            log(f"Finding optimal traversal for {contig}")
            best_score = -1
            best_delta = None
            best_path = None
            sub = []
            sub.append(ref_to_deltas[contig][0])
            for delta in sub: ## force first, remove after debug
                segs = recursive_trend_segment(delta, depth=0)
                trend = mode([seg.slope for seg in segs])
                print(f'******************************{trend}')
                find_pivots(segs, trend)
                reflect_alignment_blocks(segs, trend)
                #for seg in segs:
                #source, sink, slope = plot_trend(delta)
                source, sink, slope  = (delta.alignment_blocks[1534],delta.alignment_blocks[5649],-1)
                print(f'here {source.index}, {sink.index}')
                traversal = dijkstra_traversal(delta.alignment_blocks, source.index, sink.index, max_jump=args.max_jump, weights='euclidean')
                    
                refined_path = traversal['path_blocks']
                pdb.set_trace()
                score = score_path(refined_path, delta.reference_length)
                    #if delta.reference == '3R':
                        #pdb.set_trace()
                if score > best_score:
                    best_score = score
                    best_delta = delta
                    best_path = refined_path

                log(f"Building interval trees for {best_delta.reference}")
                domain_tree, range_tree = create_domain_range_trees(best_delta.alignment_blocks)
                log(f"Classifying SVs for {best_delta.reference}")
                new_SVs = extract_dotplot_segments_from_path(best_path, slope, best_delta.reference, domain_tree, range_tree)

                for block in best_path:
                    block.part_of_primary_synteny = True
                if args.interactive:
                    bl, tl, tr, br = get_bounding_box(delta.alignment_blocks)
                    plot_interactive_dotplot(best_delta, new_SVs, x_marker=tr[0])
            
                if repeats:
                    repeat_tree = create_repeat_tree(repeats, delta.query)
                    annotate_repeats_with_fraction(best_delta.alignment_blocks, repeat_tree, fraction_threshold=0.90)
                else:
                    repeat_tree = None

                SVs.extend(new_SVs)
                
                inversions = [b for b in best_path if b.rounded_slope * slope == -1]
                for index, inv in enumerate(inversions):
                    inv_sv = DotPlotLineSegment(
                        chrom=inv.reference,
                        reference_start=inv.reference_start,
                        reference_end=inv.reference_end,
                        query_start=inv.query_start,
                        query_end=inv.query_end,
                        sv_type='INV',
                        theta=0,
                        index=index
                    )
                    inv_sv.event_ID = f'inv_{inv.reference}_{index:04}'
                    SVs.append(inv_sv)

            if best_delta:
                if args.plot or args.interactive:
                    log(f"Plotting {best_delta.reference}")
                    if args.plot:
                        plot_static_dotplot(best_delta, new_SVs)
                    #if args.interactive:
                        #plot_interactive_dotplot(best_delta, new_SVs)

                if args.plotbw:
                    plot_bw_inset(best_delta)

        else:
            log(f"No alignments found for {contig}. Skipping...", level="WARNING")
    if args.call:
        write_vcf(SVs, output_path=args.output, sample=args.sample)
        log(f"VCF written to {args.output}")
    
    log("SVMU finished successfully.")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        log(f"Unexpected error: {str(e)}", level="ERROR")
        sys.exit(1)
