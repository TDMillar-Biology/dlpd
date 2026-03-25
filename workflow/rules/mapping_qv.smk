rule map_reads:
    input:
        reads=get_reads,
        assembly=get_assembly
    output:
        bam="results/{strain}/mapping/{strain}.bam"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping

        minimap2 -ax map-pb -t {threads} {input.assembly} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam} -O bam

        samtools index {output.bam}
        """

rule call_variants:
    input:
        bam="results/{strain}/mapping/{strain}.bam",
        assembly=get_assembly
    output:
        vcf="results/{strain}/variants/{strain}.vcf.gz"
    threads: 16
    resources:
        mem_mb=64000,
        runtime=600
    shell:
        """
        mkdir -p results/{wildcards.strain}/variants

        run_pepper_margin_deepvariant \
            -b {input.bam} \
            -f {input.assembly} \
            -o results/{wildcards.strain}/variants \
            -t {threads}
        """

rule compute_qv:
    input:
        vcf="results/{strain}/variants/{strain}.vcf.gz",
        assembly=get_assembly
    output:
        "results/{strain}/qv/qv.tsv"
    script:
        "../scripts/compute_qv.py"
