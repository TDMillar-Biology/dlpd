'''
Docstring for core.scaffold
Scaffold query contigs with reference chromosomes
Enforcing Atomicity of query space, Projection complementarity, gap exclusivity and orientation inheritence
Trevor Millar
'''
from collections import defaultdict
from intervaltree import IntervalTree
from svmu2.core.selection import select_best_projection

def atomicity_filter(alns):
    '''
    Enforce atomicity (no duplicated contigs) in query space
    '''

    query_intervals = defaultdict(IntervalTree)
    best_projections = select_best_projection(alns)
    for aln in best_projections.values(): 
        aln.compute_bounding_box()
        minx, miny = aln.bounding_box.bottom_left
        maxx, maxy = aln.bounding_box.top_right
        query_intervals[aln.query]
        domain = (minx, maxx, aln.query, aln.reference)
        arange = (miny, maxy, aln.query, aln.reference)
        print(domain)