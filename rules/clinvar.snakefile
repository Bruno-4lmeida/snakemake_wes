rule annotate_clinvar:
    input:
        vcf=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.vcf.gz",
        clinvar="/home/operator/InterOmics/workflows/interomics_snakemake/reference/clinvar.vcf.gz",
        clinvar_index="/home/operator/InterOmics/workflows/interomics_snakemake/reference/clinvar.vcf.gz.tbi",
        config="/home/operator/InterOmics/workflows/interomics_snakemake/scripts/clinvar_vcfanno.conf"
    output:
        annotated_vcf=config["output"] + "output_interomics/clinvar/{sample}/cnvs.clinvar.vcf"
    conda:
        "../envs/vcfanno.yaml"
    shell:
        """
        vcfanno {input.config} {input.vcf} > {output.annotated_vcf}
        """
