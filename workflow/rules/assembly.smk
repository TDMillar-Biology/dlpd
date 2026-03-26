rule assemble:
    input:
        reads=get_reads
    output:
        hap1="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.gfa",
        hap2="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.gfa"
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

rule gfa2fa:
    input:
        hap1="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.gfa",
        hap2="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.gfa"
    output:
        hap1fasta="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.fasta",
        hap2fasta="results/{strain}/assembly/{strain}.bp.hap2.p_ctg.fasta"
    resources:
        mem_mb=7500,
        runtime=10,
        tasks=1
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.hap1} > {output.hap1fasta}
        awk '/^S/{{print ">"$2"\\n"$3}}' {input.hap2} > {output.hap2fasta}
        """

rule mummer_for_curation:
    input:
        reference="data/reference_genomes/ISO1-r6.58_main.fasta",
        query="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.fasta"
    output:
        delta="results/{strain}/mummer/{strain}_r6_main.delta"
    resources:
        mem_mb=16000,
        runtime=120
    threads: 8
    conda:
        "../envs/mummer.yaml"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mummer

        nucmer -p results/{wildcards.strain}/mummer/{wildcards.strain}_r6_main \
            {input.reference} {input.query} -t {threads}
        """

rule curate_assembly:
    input:
        fasta="results/{strain}/assembly/{strain}.bp.hap1.p_ctg.fasta",
        breaks="config/curation/{strain}.tsv"
    output:
        fasta="results/{strain}/curated_assembly/{strain}.curated.fasta",
        log="results/{strain}/curated_assembly/{strain}.curated.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        mkdir -p results/{strain}/curated_assembly/
        breakasm \
            {input.fasta} \
            --breakfile {input.breaks} \
            --prefix results/{wildcards.strain}/curated_assembly/{wildcards.strain}
        """