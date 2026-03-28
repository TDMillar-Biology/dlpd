'''
Docstring for orchestration.parse
Orchestrate parsing of delta file and selection of main alignments
'''
from svmu2.IO.delta import parse_delta_file
from svmu2.core.selection import select_primary_alignments, select_best_projection

def run_parse(args):
    all_alns = parse_delta_file(args.delta)
    primary = select_primary_alignments(all_alns)
    #primary = select_best_projection(all_alns) ## debugging
    if not primary: 
        raise ValueError("No primary alignments found")
    
    return all_alns, primary