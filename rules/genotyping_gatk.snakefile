rule genotype_GVCFs:
    input:
        # single or list of bam files
        gvcf=config["output"] + "output_interomics/vcf/call/haplotype_caller/{sample}.g.vcf.gz",
        ref=config["refgenome"]["ref"],
        known_var=config["refgenome"]["known_var"],

    output:
        vcf=config["output"] + "output_interomics/vcf/call/haplotype_caller/{sample}.vcf.gz",
    conda:
        "../envs/gatk.yaml"
    params:
        extra="",  # optional
        java_opts="",  # optional
    threads: config["threads"]["others"],
    benchmark: "benchmarks/GATK_genotyping_{sample}_benchmark.txt"

    resources:
        mem_mb=1024,
    shell:
        """
        echo "##############################################" 
        echo "-----    Running GATK Genotyping    ----------" 
        echo "##############################################" 
               
        gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {input.gvcf} \
        -O {output.vcf} \
        --create-output-variant-index \
        -L {input.known_var}
        """
