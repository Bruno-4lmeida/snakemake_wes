rule deepvariant:
    input:
        bam= config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
        ref= config["refgenome"]["ref_dir"] + "GRCh38_numeric_genome.fasta",
    output:
        gvcf= config["output"] + "output_interomics/vcf/call/{sample}.g.vcf.gz",
        vcf= config["output"] + "output_interomics/vcf/call/{sample}.vcf.gz",
    params:
        bam_index= config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bai",
        bam_dir= config["output"] + "output_interomics/recal_apply/",
        ref_dir= config["refgenome"]["ref_dir"],
        out_dir= config["output"] + "output_interomics/vcf/call/",
        threads=16,
        model_type="WES",
    shell:
        """
        echo "##############################################" 
        echo "--------    Running DEEP VARIANT    ----------" 
        echo "##############################################"
        ./scripts/deepvariant2.sh
        """
