rule assemble:
    input:
        reads=get_reads
    output:
        primary="results/{strain}/assembly/{strain}.p_ctg.gfa",
        alternate="results/{strain}/assembly/{strain}.a_ctg.gfa"
    threads: config["threads"]
    shell:
        """
        hifiasm -t {threads} -o results/{wildcards.strain}/assembly/{wildcards.strain} {input.reads}
        """