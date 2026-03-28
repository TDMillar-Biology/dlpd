'''
Docstring for orchestration.call
Orchestrate the necessary steps for variant calling
'''

from svmu2.orchestration.synteny import run_synteny
from svmu2.core.classify import create_domain_range_trees, extract_dotplot_segments_from_path, inversion_calling
from svmu2.IO.vcf import write_vcf

## debugging
from svmu2.visualization.dotplot import render_alignment_blocks, render_aln_block, render_sv
import matplotlib.pyplot as plt
import pdb

def run_call(args):
    alns = run_synteny(args)
    SVs = []

    for ref, alignment in alns.items():
        domain_tree, range_tree = create_domain_range_trees(alignment.alignment_blocks)
        INDELS = extract_dotplot_segments_from_path(alignment.primary_synteny_blocks, alignment.slope, alignment.reference, domain_tree, range_tree)
        INVERSIONS = inversion_calling(alignment.primary_synteny_blocks, alignment.slope)
        SVs.extend(INDELS + INVERSIONS)
    
    ## Debugging
    '''
        fig, ax = render_alignment_blocks(alignment,xlabel=alignment.reference,ylabel='y')
        for b in alignment.primary_synteny_blocks:
            render_aln_block(b, ax, color = "blue")
        for sv in SVs:
            render_sv(sv, ax)
        plt.show()
        
    '''
    write_vcf(SVs, output_path=args.out, sample='SAMPLE')
    #alignment.build_primary_synteny_tree()

'''
    if args.plot_trend:
        render_trend_plot(aln, trend_result, outdir)

    if not trend_result.monotonic and args.debug:
        render_trend_debug_plot(aln, trend_result, debug_dir)


        ### The rest of this is handling the case of result.trend = "no trend" -- untold how common this is 
        # Interesting / problematic case of result.trend = "no trend"
        print("[WARN] Non-monotonic or insignificant trend detected")
        print(f"  Alignment: {aln.reference} vs {aln.query}")
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
                f"{aln.reference} vs {aln.query}\n"
                f"Trend: {result.trend}, p={result.p:.3g}"
            )
            ax.grid(True)

            outpath = os.path.join(
                debug_dir,
                f"{aln.reference}_vs_{aln.query}_trend_debug.png"
            )
            fig.savefig(outpath, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"[INFO] Debug plot written to {outpath}")
        if cumsum[0] <= cumsum[-1]:
            aln.trend = "increasing"
        else:
            aln.trend = "decreasing"
        aln.trend_checked = True
        # ----------------------------------
        # Return diagnostic payload
        # ----------------------------------
        return {
            "reference": aln.reference,
            "query": aln.query,
            "trend": result.trend,
            "p_value": result.p,
            "significant": result.h
        }
        '''