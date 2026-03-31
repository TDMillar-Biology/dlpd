rule compleasm:
    input:
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        directory("results/{strain}/compleasm")
    conda:
        "../envs/compleasm.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/compleasm/{strain}.log"
    shell:
        """
        mkdir -p {output}

        compleasm run \
            -a {input.assembly} \
            -o {output} \
            -t {threads} \
            -l diptera_odb10 \
            > {log} 2>&1
        """
