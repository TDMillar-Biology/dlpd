import argparse
import os

def main():
    ########## ARGUMENT PARSING ##########
    parser = argparse.ArgumentParser(
        prog="svmu2",
        description="Structural Variants from MUmmer and associated tools"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    ##### PLOT #####
    plot = subparsers.add_parser("plot", help="Static Dotplot of genome to genome alignment")
    plot.add_argument("--delta", '-d', required=True, help='Path to .delta alignment file')
    plot.add_argument("--out_dir", '-o', required=True, help='Out dir to write output image to')
    plot.add_argument("--strain", required=False)
    plot.add_argument("--contig", default=None)
    plot.add_argument("--reftarget", default=None, required=False, help = "plot all alignments involving reftarget")
    plot.add_argument("--qrytarget", default=None, required=False, help = "plot all alignments involving qrytarget")
    plot.add_argument("--all", action = 'store_true',default = None, required=False, help="Force plotting of all alignments")
    plot.add_argument("--popup", action = 'store_true', required = False, help="Pop up the figure as opposed to writing to file")

    ##### CALL #####
    call = subparsers.add_parser("call", help="Call SVs from genome to genome alignment. Write vcf")
    call.add_argument("--delta", '-d', required=True, help='Path to .delta alignment file')
    call.add_argument("--out", '-o', required=True, help='Path to write output vcf')
    call.add_argument('--sample', '-s', default='SAMPLE', help='Sample name for VCF output')
    call.add_argument('--write-bnds', '-wb', action='store_true', help='Also output BND breakends for each SV (default: False, niche use case, not recommended for most users)')

    ##### INTERACTIVE #####
    interactive = subparsers.add_parser("interactive", help="Interactive Dotplot of genome to genome alignment")
    interactive.add_argument("--delta", '-d', required=True, help='Path to .delta alignment file')
    interactive.add_argument("--out_dir", '-o', required=True, help='Out dir to write output html to')
    interactive.add_argument("--popup", action='store_true')

    ##### SCAFFOLD #####
    scaffold = subparsers.add_parser("scaffold", help = "Scaffold qry genome based on alignment to ref genome")
    scaffold.add_argument("--delta", '-d', required=True, help='Path to .delta alignment file')
    scaffold.add_argument('--reference','-r', help='Reference assembly for scaffolding (fasta) THIS MIGHT NOT ACTUALLY BE NEEDED =D')
    scaffold.add_argument('--query', '-q', help='Query assembly for scaffolding (fasta)')

    ##### TREND #####
    trend = subparsers.add_parser("trend", help = "Debugging function to visualize cumulative sum of weighted values of alignment blocks")
    trend.add_argument("--delta", '-d', required=True, help='Path to .delta alignment file')
    
    ##### REPEATS #####
    repeats = subparsers.add_parser("repeats", help = 'Tools for repeat annotation')
    repeats.add_argument('--repeat_masker', '-rm', help = "Path to RepeatMasker .out file")

    ### IDK WHAT TO DO WITH THESE TEMPORARILY
    '''
    parser.add_argument('--max_jump', '-mj', type=int, default=100_000, help='Maximum distance to consider for connecting blocks (see README before adjusting)')

    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')


    parser.add_argument('--reference_contigs', '-rc', action='store_true', help='Print reference contigs found in delta file')
    parser.add_argument("--canonical_contigs", help="File or comma-separated list of canonical contig names.")
    args = parser.parse_args()

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

    '''
    args = parser.parse_args()
    ########## COMMAND DISPATCH ##########

    if args.command == "plot":
        from svmu2.orchestration.plot import run_plot
        run_plot(args)

    elif args.command == "call":
        from svmu2.orchestration.call import run_call
        run_call(args)

    elif args.command == "interactive":
        from svmu2.orchestration.plot import run_plot
        run_plot(args)

    elif args.command == "scaffold":
        from svmu2.orchestration.scaffold import run_scaffold
        run_scaffold(args)

    elif args.command == "trend":
        from svmu2.orchestration.trend import run_trend
        run_trend(args)
'''
    elif args.command == "repeats":
        from svmu.repeats import run_repeats
        run_repeats(args)

    else:
        raise RuntimeError(f"Unknown command: {args.command}")

'''


if __name__ == '__main__':
    main()

