'''
Docstring for core.traversal
Produce a traversal through the graph representation of the alignment
'''    

from dataclasses import dataclass

@dataclass
class TrendResult:
    trend: str
    p_value: float
    significant: bool
    monotonic: bool

def find_source_sink_nodes(alignment, trend_result):
    """
    Determine source and sink blocks for graph traversal.
    These should be the most extreme blocks in the traversal of a monotonically increasing / decreasing plot

    returns:
        (source, sink)
    """
    # Choose corners based on trend
    if trend_result.trend == "increasing":
        start = alignment.bounding_box.bottom_left
        end   = alignment.bounding_box.top_right
    elif trend_result.trend == "decreasing":
        start = alignment.bounding_box.top_left
        end   = alignment.bounding_box.bottom_right
    else:
        raise RuntimeError(
            "Cannot determine source/sink for non-monotonic alignment"
        )

    current_min_start = None
    current_min_end = None
    min_dist_start = float("inf")
    min_dist_end = float("inf")

    for block in alignment.alignment_blocks:
        d_start = alignment.euclidean_distance(start, (block.reference_start, block.query_start))
        d_end = alignment.euclidean_distance(end, (block.reference_end, block.query_end))

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

    source = current_min_start
    sink = current_min_end

    return (source, sink)


def evaluate_alignment_trend(aln):
    import numpy as np
    import pymannkendall as mk
    aln.trace = [b.weighted_length for b in aln.alignment_blocks]
    aln.cumsum = np.cumsum(aln.trace)

    if len(aln.cumsum) < 3:
        return TrendResult(
            trend="undetermined",
            p_value=float("nan"),
            significant=False,
            monotonic=False,
        )

    result = mk.hamed_rao_modification_test(aln.cumsum, alpha=0.2)

    monotonic = result.trend in ("increasing", "decreasing")

    return TrendResult(
        trend=result.trend,
        p_value=result.p,
        significant=result.h,
        monotonic=monotonic
    )

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
    import networkx as nx
    from math import sqrt
    from sortedcontainers import SortedList

    def xdistance(a, b):
        return abs(b.reference_start - a.reference_end)

    def ydistance(a, b):
        return abs(b.query_start - a.query_end)

    def euclidean_distance(x1, y1, x2, y2):
        return sqrt((x2 - x1)**2 + (y2 - y1)**2)

    def orientation_aware_euclidean_distance(a, b):
        '''
        Accept two alignment block objects a and b
        Find the min distance between the two, allowing for reflfections
        at gany point defined on the plane (general case) including the special case of midpoint
        of the line segment
        The general case is important for complex SVs that are subsumed by inversions in reference space
        '''
        ## special case (pivot == midpoint)
        scores = []
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_start))
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.query_end))
        scores.append(euclidean_distance(a.reference_end, a.query_start, b.reference_start, b.query_start))
        scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_end, b.query_start))

        ## general case (pivot is any point on the plane)
        if getattr(a, "reflected", False):
            scores.append(euclidean_distance(a.reference_end, a.y2_reflection, b.reference_start, b.query_start)) # reflection
            #scores.append(euclidean_distance(a.reference_start, a.y1_reflection, b.reference_start, b.query_start)) # reflection
        if getattr(b, "reflected", False):
            scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_start, b.y1_reflection)) # reflection
            #scores.append(euclidean_distance(a.reference_end, a.query_end, b.reference_end, b.y2_reflection)) # reflection
        return min(scores)
    
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
        'euclidean': orientation_aware_euclidean_distance,
        'X': xdistance,
        'Y': ydistance
    }.get(weights)

    if weight_fn is None:
        raise ValueError(f"Invalid weight designation for Dijkstra traversal: {weights}")

    # Leave it to the user to define a good max distance -- leaving this structure for expansion
    for factor in [1]:
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

def compute_primary_synteny(alignment, max_jump=100_000, weights="euclidean"):
    if alignment.source is None or alignment.sink is None:
        raise RuntimeError(
            f"compute_primary_synteny called on non-traversable alignment "
            f"{alignment.reference} vs {alignment.query}"
        )

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