rule filter_var:
    input:
        list(filter(None, [
            config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_snps.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_snps.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_snps.vcf.gz" if config.get("use_deepvariant", False) else None,
            config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_indels.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_indels.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_indels.vcf.gz" if config.get("use_deepvariant", False) else None
        ])),
    output:
        list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_snps.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_snps.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_snps.vcf.gz" if config.get("use_deepvariant", False) else None,
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_indels.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_indels.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_indels.vcf.gz" if config.get("use_deepvariant", False) else None
        ])),
    params:
        ref=config["refgenome"]["ref"],
    conda:
        "../envs/gatk.yaml",
    shell:
        """
        echo "##############################################" 
        echo "-----   GATK Variant Filter    --------------" 
        echo "##############################################"

        for vcf in {input}; do
            if [ -f $vcf ]; then
                echo "Processing $vcf"
                output_file=$(echo $vcf | sed 's/selected/filtered/' | sed 's/raw/filtered/')

                gatk VariantFiltration -R {params.ref} -V $vcf \
                    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                    --filter-name "variant_filter" -G-filter "GQ < 20.0" -G-filter-name "lowGQ" \
                    -O $output_file &
            fi
        done

        wait
        """
