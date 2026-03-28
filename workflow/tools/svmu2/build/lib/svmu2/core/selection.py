'''
Docstring for core.selection
select the main alignments found in delta file
'''

from collections import defaultdict

def union_intervals(intervals):
    '''
    Docstring for union_intervals
    
    let Qi = some interval (x1i, x2i) on a line segment bound by I = (X1, X2)
    s.t. X1 <= x1i < x2i <= X2
    and Q = set(Q1 ... Qn)
    Then there should exist UnionQ, defined by Q1 U Q2 U ... Qi
    Accept Q, return UnionQ
    '''
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = []

    cur_start, cur_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end

    merged.append((cur_start, cur_end))
    return merged

def complement(intervals, X1, X2):
    '''
    let Qi = some interval (x1i, x2i) on a line segment bound by I = (X1, X2)
    s.t. X1 <= x1i < x2i <= X2
    and Q = set(Q1 ... Qn)
    Q must be sorted by start position and there must be no overlaps such that there is no
    Qi U Qk != 0 for all i, k where i!=k
    
    really just run this after you compute the union of the intervals and you'll be safe
    '''
    if not intervals:
        return [(X1, X2)]

    comp = []
    cur = X1

    for start, end in intervals:
        if start > cur:
            comp.append((cur, start))
        cur = end

    if cur < X2:
        comp.append((cur, X2))

    return comp

def interval_length(intervals):
    return sum(end - start for start, end in intervals)

def select_primary_alignments(alignments):
    '''
    1-1 mapping of ref contig to best qry contig
    accept iterable of alignments, 
    return dictionary of refchrom:qrychrom for the best alns
    '''
    by_ref = defaultdict(list)
    for aln in alignments:
        by_ref[aln.reference].append(aln)

    primary = {}

    for ref, alns in by_ref.items():
        best = None
        best_cov = -1

        for aln in alns:
            intervals = []
            for block in aln.alignment_blocks: # I realize this is very defensive for ref space
                start = min(block.reference_start, block.reference_end)
                end = max(block.reference_start, block.reference_end)
                intervals.append((start, end))

            union = union_intervals(intervals)
            cov = interval_length(union)

            if cov > best_cov:
                best_cov = cov
                best = aln

        primary[ref] = best

    return primary

def select_best_projection(alignments):
    """
    1-1 mapping of query contig to best reference contig.
    Best defined as maximal query-space coverage.
    """
    by_qry = defaultdict(list)
    for aln in alignments:
        by_qry[aln.query].append(aln)

    primary = {}

    for qry, alns in by_qry.items():
        best = None
        best_score = -1

        for aln in alns:
            intervals = []
            for block in aln.alignment_blocks:
                start = min(block.query_start, block.query_end)
                end = max(block.query_start, block.query_end)
                intervals.append((start, end))

            union = union_intervals(intervals)
            cov = interval_length(union)

            score = cov #/ aln.query_length
            aln.score = score
            if score > best_score:
                best_score = score
                best = aln

        primary[qry] = best

    return primary