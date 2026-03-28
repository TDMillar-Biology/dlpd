'''
core.classify

Operations to classify SV types according to geometry of models.line_segment objects
'''

from svmu2.models.line_segment import DotPlotLineSegment
import pdb

def is_small(x, threshold = 20): #20 nucleotides is what mummer is capable of detecting per publication
    return abs(x) < threshold

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

def create_domain_range_trees(alignment_blocks, min_overlap_buffer=25):
    '''
    Accept Aligment blocks as iterable
    ensures start < end as interval tree requires
    returns tuple of objects being (domain_tree, range_tree)
    This implementation of an interval tree allows me to search for overlapping segments in log n time
    ''' 
    from intervaltree import IntervalTree
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

def find_theta(obj1, obj2):
    ''' 
    Calculate the angle theta (in degrees) in a polar coordinate plane 
    that defines the direction of the line segment from obj1 to obj2.
    '''
    import math
    # Calculate change in x and y
    delta_x = obj2.reference_start - obj1.reference_end
    delta_y = obj2.query_start - obj1.query_end

    # Use atan2 for proper angle computation across all quadrants
    angle_rad = math.atan2(delta_y, delta_x) # div by zero protection returns (-pi,pi) 
    angle_deg = math.degrees(angle_rad)
    
    return angle_deg

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

def inversion_calling(syntenic_path, slope):
    inversion_svs = []
    inversions = [b for b in syntenic_path if b.rounded_slope * slope == -1]
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

        inversion_svs.append(inv_sv)
    return inversion_svs