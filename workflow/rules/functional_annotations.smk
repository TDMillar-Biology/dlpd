rule prepare_functional_gff:
    input:
        gff="data/reference_genomes/dmel-all-r6.58.gff"
    output:
        gff="data/reference_genomes/dmel-r6.58.functional.gff"
    resources:
        mem_mb=1000,
        runtime=15,
        ntasks=1
    container:
        "workflow/containers/images/genomic_features.sif"
    log:
        "logs/variants/prepare_functional_gff.log"
    shell:
        r"""
        mkdir -p logs/variants
        exec 2> {log:q}

        awk -F '\t' '
            BEGIN {{ print "##gff-version 3" }}
            /^##FASTA/ {{ exit }}
            /^#/ {{ next }}
            NF == 9 {{
                text = tolower($2 "\t" $3 "\t" $9)
                if (text ~ /chip_seq|redfly|promoter/ || $3 == "TF_binding_site" || $3 == "regulatory_region" || $3 == "exon")
                    print
            }}
        ' {input.gff:q} > {output.gff:q}
        """

rule subset_sv_euchromatin:
    input:
        vcf="results/{strain}/variants/{strain}.svmu2.vcf",
        bed="config/euchromatin_boundaries.bed"
    output:
        vcf="results/{strain}/variants/{strain}.euchromatic.svmu2.vcf.gz",
        index="results/{strain}/variants/{strain}.euchromatic.svmu2.vcf.gz.csi"
    resources:
        mem_mb=4000,
        runtime=15,
        ntasks=1
    container:
        "workflow/containers/images/genomic_features.sif"
    log:
        "logs/variants/{strain}.subset_euchromatin.log"
    shell:
        r"""
        mkdir -p results/{wildcards.strain}/variants logs/variants
        exec > {log:q} 2>&1

        bcftools view \
            --targets-file {input.bed:q} \
            --targets-overlap 1 \
            --output-type u \
            {input.vcf:q} | \
        bcftools sort \
            --max-mem 1G \
            --output-type z \
            --output {output.vcf:q}

        bcftools index --force --csi {output.vcf:q}
        """

rule intersect_sv_functional_annotations:
    input:
        vcf="results/{strain}/variants/{strain}.euchromatic.svmu2.vcf.gz",
        index="results/{strain}/variants/{strain}.euchromatic.svmu2.vcf.gz.csi",
        gff="data/reference_genomes/dmel-r6.58.functional.gff"
    output:
        overlaps="results/{strain}/variants/{strain}.sv_gff_overlaps.tsv"
    resources:
        mem_mb=16000,
        runtime=60,
        ntasks=1
    container:
        "workflow/containers/images/genomic_features.sif"
    log:
        "logs/variants/{strain}.functional_intersect.log"
    shell:
        r"""
        mkdir -p results/{wildcards.strain}/variants logs/variants
        exec 2> {log:q}

        bedtools intersect \
            -a {input.vcf:q} \
            -b {input.gff:q} \
            -wa -wb \
            > {output.overlaps:q}
        """
