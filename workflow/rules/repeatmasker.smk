# Include this file from the workflow Snakefile.
# Input chromosome names must be 2L, 2R, 3L, 3R, 4, and X in both FASTA and AGP.

rule repeatmask_scaffolded:
    input:
        fasta="results/{strain}/scaffold/{strain}.scaffolded.fasta",
        library="data/repeat_library/LarracuenteLab.Repeat.library.specieslib_mod2_Mel_042623.fasta"
    output:
        rm_out="results/{strain}/repeats/{strain}.scaffolded.fasta.out",
        normalized="results/{strain}/repeats/{strain}.rm.tsv"
    params:
        outdir="results/{strain}/repeats"
    threads: 
    resources:
        mem_mb=32000,
        runtime=720,
        ntasks=1
    container:
        "workflow/containers/images/repeatmasker.sif"
    log:
        "logs/repeats/{strain}.repeatmasker.log"
    shell:
        r"""
        mkdir -p {params.outdir:q} logs/repeats
        exec > {log:q} 2>&1

        # Each parallel RMBlast batch uses four cores.
        # Snakemake may reduce threads to the available core allocation.
        if [ {threads} -lt 4 ]; then
            echo "RepeatMasker with RMBlast needs at least 4 allocated cores." >&2
            exit 1
        fi
        parallel_jobs=$(( {threads} / 4 ))

        RepeatMasker -v
        RepeatMasker \
            -engine rmblast \
            -famdb_dir "" \
            -lib {input.library:q} \
            -pa "$parallel_jobs" \
            -dir {params.outdir:q} \
            {input.fasta:q}

        rmtools normalize \
            --rm-out {output.rm_out:q} \
            --out {output.normalized:q} \
            --strain {wildcards.strain:q}
        """


rule plot_repeat_panels:
    input:
        rm="results/{strain}/repeats/{strain}.rm.tsv",
        agp="results/{strain}/scaffold/{strain}.scaffolded.agp"
    output:
        plots=expand(
            "results/{{strain}}/repeats/plots/{chrom}.panel.pdf",
            chrom=["2L", "2R", "3L", "3R", "4", "X"]
        )
    params:
        outdir="results/{strain}/repeats/plots",
        taxonomy="class",
        rm_bin=50000
    threads: 1
    resources:
        mem_mb=8000,
        runtime=120,
        ntasks=1
    container:
        "workflow/containers/images/repeatmasker.sif"
    log:
        "logs/repeats/{strain}.panels.log"
    shell:
        r"""
        mkdir -p {params.outdir:q} logs/repeats
        exec > {log:q} 2>&1

        # Both tracks are in scaffold coordinates. AGP supplies the full extent,
        # including sequence beyond the last repeat annotation.
        for chrom in 2L 2R 3L 3R 4 X; do
            rmtools panel \
                --rm {input.rm:q} \
                --agp {input.agp:q} \
                --region "$chrom" \
                --taxonomy {params.taxonomy:q} \
                --rm-bin {params.rm_bin} \
                --out {params.outdir:q}/"$chrom.panel.pdf"
        done
        """
