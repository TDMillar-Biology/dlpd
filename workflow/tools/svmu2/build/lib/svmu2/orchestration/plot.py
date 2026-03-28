'''
Orchestrate the necessary steps for the plot command
'''

from svmu2.models.line import LinePrimitive, build_alignment_primitives
from svmu2.orchestration.parse import run_parse
import os
from pathlib import Path


#debugging
import matplotlib.pyplot as plt
def run_plot(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    all_alns, primary = run_parse(args)

    targets = primary.values()
    if getattr(args, "reftarget", None):
        targets = [aln for aln in all_alns if aln.reference == args.reftarget]
        if not targets:
            raise ValueError(f"{args.reftarget} not a reference sequence in {args.delta}")
    if getattr(args, "qrytarget", None):
        targets = [aln for aln in all_alns if aln.query == args.qrytarget]
        if not targets:
            raise ValueError(f"{args.qrytarget} not a query sequence in {args.delta}")
    if getattr(args, "all", None):
        targets = all_alns

    for aln in targets:
        primitives = build_alignment_primitives(aln)

        xlabel = aln.reference
        ylabel = aln.query
        
        if args.command == 'interactive':
            from svmu2.visualization.renderers import render_plotly
            outpath = out_dir / f"{aln.reference}_{aln.query}.html"
            title = f"{aln.reference}_{aln.query}"
            fig = render_plotly(primitives, title=title)
            
            if getattr(args, "popup", None):
                fig.show()
            else:
                fig.write_html(outpath)
        else:
            from svmu2.visualization.renderers import render_matplotlib
            fig, ax = render_matplotlib(primitives, xlabel, ylabel)
            outpath = out_dir / f"{aln.reference}.{aln.query}.pdf"
            if getattr(args, "popup", None):
                plt.show()
            else:
                fig.savefig(outpath)
        plt.close()
