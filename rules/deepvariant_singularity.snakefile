rule deepvariant:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
        ref=config["refgenome"]["ref"] 
    output:
        vcf=config["output"] + "output_interomics/vcf/call/deep_variant/{sample}.vcf.gz",
    threads: config["threads"]["deepvariant"],
    params:
        bam_dir=config["output"] + "output_interomics/recal_apply/",
        ref_dir=config["refgenome"]["ref_dir"],
        bed=config["refgenome"]["bed"],
        out_dir=config["output"] + "output_interomics/vcf/call/deep_variant/",
        threads=config["threads"]["deepvariant"],
        model_type=config["model_type"],
        deep_sif=config["deepvariant_sif_path"],
    benchmark: config["output"] + "output_interomics/benchmarks/caller/deep_variant/deepvariant_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################"
        echo "--------    Running DEEP VARIANT    ----------"
        echo "##############################################"

        singularity exec \
            --bind {params.bam_dir}:/input \
            --bind {params.ref_dir}:/reference \
            --bind {params.out_dir}:/output \
            {params.deep_sif} \
            /opt/deepvariant/bin/run_deepvariant \
            --num_shards={params.threads} \
            --model_type={params.model_type} \
            --ref={input.ref} \
            --reads=/input/{wildcards.sample}_markduplicates_recal.bam \
            --output_vcf=/output/{wildcards.sample}.vcf.gz \
            --regions={params.bed} 
        """
