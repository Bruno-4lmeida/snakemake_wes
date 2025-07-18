rule funcotator:
    input:
        vcf=[
            (config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz" if config["use_gatk"] == True else []),
            (config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz" if config["use_strelka"] == True else []),
            (config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz" if config["use_deepvariant"] == True else [])
        ],
        ref=config["refgenome"]["ref"]
    output:
        [
            (config["output"] + "output_interomics/funcotator/haplotype_caller/{sample}_funcotator.vcf.gz" if config["use_gatk"] == True else []),
            (config["output"] + "output_interomics/funcotator/strelka/{sample}_funcotator.vcf.gz" if config["use_strelka"] == True else []),
            (config["output"] + "output_interomics/funcotator/deep_variant/{sample}_funcotator.vcf.gz" if config["use_deepvariant"] == True else []),
        ]
    conda:
        "../envs/gatk.yaml"
    params:
        data_source=config["refgenome"]["data_source"]
    threads: config["threads"]["others"]
    resources:
        mem_mb=1024
    benchmark:
        [
            (config["output"] + "output_interomics/benchmarks/annotated/haplotype_caller/funcotator_{sample}_benchmark.txt" if config["use_gatk"] == True else []),
            (config["output"] + "output_interomics/benchmarks/annotated/strelka/funcotator_{sample}_benchmark.txt" if config["use_strelka"] == True else []),
            (config["output"] + "output_interomics/benchmarks/annotated/deep_variant/funcotator_{sample}_benchmark.txt" if config["use_deepvariant"] == True else [])
        ]
    shell:
        """
        echo "##############################################" 
        echo "--------   Running GATK Funcotator  ---------" 
        echo "##############################################"
        gatk Funcotator \
        -R {input.ref} \
        -V {input.vcf} \
        -O {output} \
        --output-file-format VCF \
        --data-sources-path {params.data_source}/ \
        --ref-version hg38
        """
