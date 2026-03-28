'''
Docstring for orchestration.synteny
Operate on alignment graph space to determine elements of primary synteny
'''
from svmu2.orchestration.parse import run_parse
from svmu2.core.traversal import evaluate_alignment_trend, find_source_sink_nodes, dijkstra_traversal
import pdb
def run_synteny(args):
    all_alns, primary = run_parse(args)
    for ref, alignment in primary.items(): 

        # Exhaustive homologous blocks
        #build_exhaustive_tree(alignment)

        # Trend + geometry
        trend_result = evaluate_alignment_trend(alignment)

        if trend_result.monotonic:
            if trend_result.trend == 'increasing':
                alignment.slope = 1
            if trend_result.trend == 'decreasing':
                alignment.slope = -1

            alignment.compute_bounding_box()
            source, sink = find_source_sink_nodes(alignment, trend_result)
            source_idx = alignment.alignment_blocks.index(source)
            sink_idx   = alignment.alignment_blocks.index(sink)

        # Primary synteny traversal
            traversal = dijkstra_traversal(
                alignment.alignment_blocks,
                source_idx,
                sink_idx,
                max_jump=1_000_000,
                weights='euclidean',
            )
            syntenic_path = traversal["path_blocks"]

            alignment.set_primary_synteny_blocks(syntenic_path)
    return primary