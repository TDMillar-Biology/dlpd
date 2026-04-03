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

        compleasm.py run \
            -a {input.assembly} \
            -o {output} \
            -t {threads} \
            -l diptera \
            --odb odb10 \
            -L data/ortholog_databases \
            > {log} 2>&1
        """

rule busco:
    input:
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        directory("results/{strain}/busco")
    conda:
        "../envs/busco.yaml"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=120
    log:
        "logs/busco/{strain}.log"
    params:
        lineage=config["busco"]["lineage"],
        db_dir=config["busco"]["db_dir"]
    shell:
        """
        busco \
            -i {input.assembly} \
            -o busco \
            --out_path results/{wildcards.strain} \
            -l {params.lineage} \
            --download_path {params.db_dir} \
            -m genome \
            -c {threads} \
            --offline \
            -f \
            > {log} 2>&1
        """

rule compute_qv:
    input:
        vcf="results/{strain}/variants/{strain}_pmdv.vcf.gz",
        assembly="results/{strain}/assembly/{strain}.bp.p_ctg.fasta"
    output:
        dir=directory("results/{strain}/qc"),
        tsv="results/{strain}/qc/{strain}_qc.tsv"
    conda:
        "../envs/assembly_qc.yaml"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1
    log:
        "logs/qv/{strain}.log"
    shell:
        """
        mkdir -p "{output.dir}"

        python3 workflow/scripts/qv_from_vcf.py \
            --vcf {input.vcf} \
            --strain {wildcards.strain} \
            --fasta {input.assembly} \
            --output {output.tsv} \
            --min_gq 20 \
            > {log}
        """

