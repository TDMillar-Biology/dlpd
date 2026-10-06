rule align_iso1_to_assembly:
    input:
        target_assembly="results/{strain}/curated_assembly/{strain}.curated.fasta",
        query_iso1=config["references"]["ISO1"]["full"]
    output:
        paf="results/{strain}/mapping/ISO1_to_{strain}.raw.paf"
    threads: 8
    resources:
        mem_mb=16000,
        runtime=60,
        ntasks=1,
        slurm_partition="medium"
    container:
        "workflow/containers/images/mapping_qc.sif"
    log:
        "logs/mapping/ISO1_to_{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/mapping logs/mapping

        # Direct PAF emission with --cs for paftools.js compatibility
        minimap2 -cx asm5 --cs --secondary=no -t {threads} \
            {input.target_assembly} {input.query_iso1} \
            > {output.paf} 2> {log}
        """

rule svmu_filter_main_diagonal:
    input:
        paf="results/{strain}/mapping/ISO1_to_{strain}.raw.paf"
    output:
        filtered_paf=temp("results/{strain}/mapping/ISO1_to_{strain}.filtered.paf"),
        sorted_paf="results/{strain}/mapping/ISO1_to_{strain}.filtered.sorted.paf"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1,
        slurm_partition="short"
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/mapping/svmu_filter_{strain}.log"
    shell:
        """
        # 1. Filter PAF to the main diagonal using native PAF format
        svmu2 filter \
            -a {input.paf} \
            --format paf \
            -o {output.filtered_paf} \
            > {log} 2>&1

        # 2. Sort directly by Target Name (k6) and Target Start (k8) for paftools.js
        sort -k6,6 -k8,8n {output.filtered_paf} > {output.sorted_paf} 2>> {log}
        """

rule paftools_call_micro_variants:
    """
    Calls small indels and SNPs from the purified, sorted syntenic PAF.
    """
    input:
        target_assembly="results/{strain}/curated_assembly/{strain}.curated.fasta",
        paf="results/{strain}/mapping/ISO1_to_{strain}.filtered.sorted.paf"
    output:
        vcf="results/{strain}/mappable_variants/ISO1_to_{strain}.paftools_micro.vcf"
    threads: 8
    resources:
        mem_mb=8000,
        runtime=30,
        ntasks=1,
        slurm_partition="short"
    container:
        "workflow/containers/images/mapping_qc.sif"
    log:
        "logs/variants/paftools_call_{strain}.log"
    shell:
        """
        mkdir -p results/{wildcards.strain}/variants logs/variants
        paftools.js call -f {input.target_assembly} {input.paf} > {output.vcf} 2> {log}
        """

rule svmu_call:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_scaffolded.delta"
    output:
        vcf="results/{strain}/variants/{strain}.svmu2.vcf",
        bedpe="results/{strain}/variants/{strain}.svmu2.bedpe",
        plots=directory("results/{strain}/variants/svmu2_call_plots")
    resources:
        mem_mb=128000,
        runtime=60,
        ntasks=1,
        slurm_partition="short"
    log:
        "logs/variants/svmu2_call_{strain}.log"
    container:
        "workflow/containers/images/svmu2.sif"
    shell:
        """
        mkdir -p results/{wildcards.strain}/variants logs/variants

        svmu2 call \
            --alignment {input.delta} \
            --format delta \
            --out {output.vcf} \
            --sample {wildcards.strain} \
            --write-bedpe \
            --plot-svs \
            > {log} 2>&1
        """

rule summarize_sv_bedpe:
    input:
        bedpe="results/{strain}/variants/{strain}.svmu2.bedpe",
        script="workflow/scripts/summarize_svs.py"
    output:
        chromosome="results/{strain}/variants/{strain}.chromosome_summary.tsv",
        sample="results/{strain}/variants/{strain}.sample_summary.tsv"
    resources:
        mem_mb=4000,
        runtime=15,
        ntasks=1
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        python {input.script} \
            --bedpe {input.bedpe} \
            --sample {wildcards.strain} \
            --out-prefix results/{wildcards.strain}/variants/{wildcards.strain}
        """

rule svmu_plot:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_scaffolded.delta",
        ref_bed="config/euchromatin_boundaries.bed",
        qry_bed="results/{strain}/euchromatin/{strain}.euchromatin.bed"
    output:
        plots=directory("results/{strain}/variants/svmu2_synteny_plots")
    resources:
        mem_mb=32000,
        runtime=60,
        ntasks=1,
        slurm_partition="short"
    log:
        "logs/variants/svmu2_plot_{strain}.log"
    container:
        "workflow/containers/images/svmu2.sif"
    shell:
        r"""
        mkdir -p {output.plots:q} logs/variants

        svmu2 plot \
            --alignment {input.delta:q} \
            --format delta \
            --out_dir {output.plots:q} \
            --img-format pdf \
            --synteny \
            --ref-bed {input.ref_bed:q} \
            --qry-bed {input.qry_bed:q} \
            > {log:q} 2>&1
        """

rule plot_structural_geometry:
    input:
        geometry="results/{strain}/variants/{strain}.sv_geometry.tsv",
        script="workflow/scripts/plot_structural_geometry.py"
    output:
        span="results/{strain}/variants/{strain}.span_geometry.png",
        theta="results/{strain}/variants/{strain}.theta_distribution.png"
    resources:
        mem_mb=4000,
        runtime=15,
        ntasks=1
    container:
        "workflow/containers/images/python_mummer.sif"
    shell:
        """
        python {input.script} \
            --input {input.geometry} \
            --out-prefix results/{wildcards.strain}/variants/{wildcards.strain}
        """

rule project_euchromatin:
    input:
        delta="results/{strain}/mummer/{strain}_r6_main_scaffolded.delta",
        bed="config/euchromatin_boundaries.bed",
        script="workflow/scripts/project_euchromatin.py"
    output:
        bed="results/{strain}/euchromatin/{strain}.euchromatin.bed",
        boundaries="results/{strain}/euchromatin/{strain}.boundaries.tsv"
    params:
        prefix="results/{strain}/euchromatin/{strain}"
    threads: 1
    resources:
        mem_mb=50000,
        runtime=60,
        ntasks=1
    container:
        "workflow/containers/images/svmu2.sif"
    log:
        "logs/euchromatin/{strain}.log"
    shell:
        r"""
        mkdir -p results/{wildcards.strain}/euchromatin logs/euchromatin

        python {input.script:q} \
            --delta {input.delta:q} \
            --bed {input.bed:q} \
            --out-prefix {params.prefix:q} \
            > {log:q} 2>&1
        """