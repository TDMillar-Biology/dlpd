rule map_reads:
    input:
        reads="data/{strain}.fq.gz",
        assembly="data/{strain}.fa"
    output:
        bam="results/{strain}/mapping/aln.bam"
    threads: config["threads"]
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        minimap2 -ax map-pb -t {threads} {input.assembly} {input.reads} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule compute_depth:
    input:
        bam="results/{strain}/mapping/aln.bam"
    output:
        depth="results/{strain}/depth/depth.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        samtools depth -a {input.bam} > {output.depth}
        """

rule compute_qv:
    input:
        depth="results/{strain}/depth/depth.txt"
    output:
        "results/{strain}/qv/qv.tsv"
    script:
        "../scripts/compute_qv.py"