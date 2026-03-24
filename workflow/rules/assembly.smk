rule assemble:
    input:
        reads=get_reads
    output:
        primary="results/{strain}/assembly/{strain}.p_ctg.gfa",
        alternate="results/{strain}/assembly/{strain}.a_ctg.gfa"
    threads: 16
    conda:
        "../envs/assembly.yaml"
    resources:
        mem_mb=64000,
        runtime=600,
        tasks=1
    shell:
        """
        mkdir -p results/{wildcards.strain}/assembly

        hifiasm -t {threads} \
            -o results/{wildcards.strain}/assembly/{wildcards.strain} \
            {input.reads} \
            > results/{wildcards.strain}/assembly/hifiasm.log 2>&1
        """

