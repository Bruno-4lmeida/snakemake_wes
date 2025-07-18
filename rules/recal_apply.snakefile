rule recal_apply:
    input:
        bam=config["output"] + "output_interomics/marked/{sample}_markduplicates.bam",
        ref=config["refgenome"]["ref"],
        recal_table=config["output"] + "output_interomics/recal_table/{sample}_recal.grp",
    output:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
    conda:
        "../envs/gatk.yaml"
    threads: config["threads"]["recal_apply"],
    benchmark: config["output"] + "output_interomics/benchmarks/recal_apply/apply_table_{sample}_benchmark.txt",
    shell:
        """
        echo "##############################################" 
        echo "-------    Running GATK AplyBQSR    ---------" 
        echo "##############################################"

        gatk ApplyBQSR \
        -I {input.bam} \
        -R {input.ref} \
        --bqsr-recal-file {input.recal_table} \
        -O {output.bam}
        """