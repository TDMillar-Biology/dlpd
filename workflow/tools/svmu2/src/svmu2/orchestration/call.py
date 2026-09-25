'''
Docstring for orchestration.call
Orchestrate the necessary steps for variant calling
'''

from pathlib import Path

from svmu2.orchestration.synteny import run_synteny
from svmu2.core.classify import (
    create_domain_range_trees,
    extract_collinear_gap_segments_from_path,
    call_all_inversions,
)
from svmu2.IO.vcf import write_vcf

class IntraChainVariant:
    """Lightweight object to perfectly mimic a DotPlotLineSegment for write_vcf"""
    def __init__(self, sv_dict):
        self.chrom = sv_dict.get("chrom")
        self.sv_type = sv_dict.get("svtype")

        self.reference_start = sv_dict.get("pos")
        self.reference_end = sv_dict.get("end")

        # write_vcf uses (query_end - query_start) to calculate INS length.
        # We can mock this by setting start to 0 and end to the actual length.
        svlen = sv_dict.get("svlen")
        self.query_start = 0
        self.query_end = svlen

        # Default structural attributes expected by write_vcf logic
        self.theta = 0
        self.range_partners = None
        self.domain_partners = None
        self.event_ID = f"intra_{self.reference_start}_{self.sv_type}"


def write_bedpe6(SVs, output_path):
    count = 0

    with open(output_path, "w") as out_f:
        for sv in SVs:
            # Intra-chain variants currently contain synthetic query
            # coordinates used only for VCF writing, so they cannot yet
            # be represented correctly as BEDPE6.
            if isinstance(sv, IntraChainVariant):
                continue

            ref_start, ref_end = sorted(
                (sv.reference_start, sv.reference_end)
            )
            query_start, query_end = sorted(
                (sv.query_start, sv.query_end)
            )

            out_f.write(
                f"{sv.reference_chrom}\t"
                f"{ref_start}\t"
                f"{ref_end}\t"
                f"{sv.query_chrom}\t"
                f"{query_start}\t"
                f"{query_end}\n"
            )

            count += 1

    print(
        f"BEDPE6 export complete. "
        f"Wrote {count:,} SVs -> {output_path}"
    )

def call(args):
    alns = run_synteny(args)
    SVs = []

    if getattr(args, "plot_svs", False):
        from svmu2.visualization.renderers import plot_interactive_sv_calls

        plot_dir = Path(args.out).parent / "svmu2_call_plots"
        plot_dir.mkdir(parents=True, exist_ok=True)

    for _, alignment in alns.items():
        domain_tree, range_tree = create_domain_range_trees(
            alignment.alignment_blocks
        )

        # 1. Call inter-chain collinear gaps (Large INDELs)
        INDELS = extract_collinear_gap_segments_from_path(
            alignment.primary_synteny_blocks,
            alignment.reference,
            domain_tree,
            range_tree,
            write_bnds=args.write_bnds,
        )

        # 2. Call inter-chain inversions
        INVERSIONS = call_all_inversions(
            alignment.final_path_segments,
            alignment.primary_synteny_blocks,
            alignment.slope,
        )

        alignment_SVs = INDELS + INVERSIONS

        for sv in alignment_SVs:
            sv.reference_chrom = alignment.reference
            sv.query_chrom = alignment.query

        SVs.extend(alignment_SVs)

        # Optional interactive debugging plot
        if getattr(args, "plot_svs", False):
            plot_path = (
                plot_dir /
                f"{alignment.reference}.{alignment.query}.html"
            )

            fig = plot_interactive_sv_calls(
                alignment,
                alignment_SVs,
                plot_path,
                auto_open=False,
            )

            del fig

        # 3. Call intra-chain micro-indels
        if getattr(args, "include_intra", False):
            intra_variants = []

            if alignment.primary_synteny_blocks:
                for block in alignment.primary_synteny_blocks:
                    block_indels = block.call_intra_chain_indels()

                    for sv_dict in block_indels:
                        intra_variants.append(
                            IntraChainVariant(sv_dict)
                        )

            SVs.extend(intra_variants)

    return alns, SVs


def run_call(args):
    alns, SVs = call(args)

    write_vcf(
        SVs,
        output_path=args.out,
        sample=args.sample,
    )

    if args.write_bedpe:
        bedpe_path = str(Path(args.out).with_suffix(".bedpe"))
        write_bedpe6(SVs, bedpe_path)

