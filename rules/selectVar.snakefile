rule select_var:
    input:
        list(filter(None, [
            config["output"] + "output_interomics/vcf/call/haplotype_caller/{sample}.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/call/strelka/{sample}.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/call/deep_variant/{sample}.vcf.gz" if config.get("use_deepvariant", False) else None
        ])),
    output:
        list(filter(None, [
            config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_snps.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_indels.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_indels.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_snps.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_snps.vcf.gz" if config.get("use_deepvariant", False) else None,
            config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_indels.vcf.gz" if config.get("use_deepvariant", False) else None
        ])),
    params:
        ref=config["refgenome"]["ref"],
    conda:
        "../envs/gatk.yaml",
    threads:5
    shell:
        """
        echo "##############################################" 
        echo "-----   GATK Select and Filter Variants ------" 
        echo "##############################################"

        for vcf in {input}; do
            if [ -f $vcf ]; then
                echo "Processing $vcf"
                output_snps=$(echo $vcf | sed 's/call/selected/' | sed 's/.vcf.gz/_raw_snps.vcf.gz/')
                output_indels=$(echo $vcf | sed 's/call/selected/' | sed 's/.vcf.gz/_raw_indels.vcf.gz/')

                gatk SelectVariants -R {params.ref} -V $vcf --select-type-to-include SNP -O $output_snps &
                gatk SelectVariants -R {params.ref} -V $vcf --select-type-to-include INDEL -O $output_indels &
            fi
        done

        # Esperar todos os processos paralelos terminarem
        wait
        """
