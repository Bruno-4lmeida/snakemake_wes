rule haplotype_caller:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
        ref=config["refgenome"]["ref"],
        dict=config["refgenome"]["dict"],
        known_var=config["refgenome"]["known_var"],
        bed=config["refgenome"]["bed"],
    output:
        vcf=config["output"] + "output_interomics/vcf/call/haplotype_caller/{sample}.vcf.gz",
    conda:
        "../envs/gatk.yaml"
    params:
        extra="",  # optional
        java_opts="",  # optional
    threads: config["threads"]["others"],
    resources:
        mem_mb=1024,
    shell:
        """
        echo "##############################################" 
        echo "-----  Running GATK Haplotype Caller   -------" 
        echo "##############################################"  

        gatk HaplotypeCaller  \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.vcf} \
        --dbsnp {input.known_var} \
        -L {input.bed}
        
        """

