rule map_reads:
    input:
        reads=get_reads,
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        bam="results/{strain}/mapping/{strain}.bam",
        bai="results/{strain}/mapping/{strain}.bam.bai"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600,
        tasks=1
    conda:
        "../envs/mapping.yaml"
    log:
        "logs/mapping/{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping

        minimap2 -ax map-pb -t {threads} {input.assembly} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam} -O bam

        samtools index {output.bam} {output.bai} \
        > {log} 2>&1
        """

rule call_variants:
    input:
        bam="results/{strain}/mapping/{strain}.bam",
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        outdir=directory("results/{strain}/variants")
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600
    log:
        "logs/variants/{strain}.log"
    singularity:
        "../containers/pmdv.sif"
    shell:
        """
        mkdir -p {output.outdir}

        run_pepper_margin_deepvariant \
            -b {input.bam} \
            -f {input.assembly} \
            -o {output.outdir} \
            -t {threads} \
            > {log} 2>&1
        """
rule compute_qv:
    input:
        vcf="results/{strain}/variants/{strain}.vcf.gz",
        assembly=get_assembly
    output:
        "results/{strain}/qv/qv.tsv"
    script:
        "../scripts/compute_qv.py"
