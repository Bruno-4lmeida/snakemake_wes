rule concat_bcftools:
    input:
        snp=list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_snps.vcf.gz" if config["use_gatk"] else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_snps.vcf.gz" if config["use_strelka"] else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_snps.vcf.gz" if config["use_deepvariant"] else None,
        ])),
        indel=list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_indels.vcf.gz" if config["use_gatk"] else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_indels.vcf.gz" if config["use_strelka"] else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_indels.vcf.gz" if config["use_deepvariant"] else None,
        ])),
    output:
        out=list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered.vcf.gz" if config["use_gatk"] else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered.vcf.gz" if config["use_strelka"] else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered.vcf.gz" if config["use_deepvariant"] else None,
        ])),
        out_final=list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz" if config["use_gatk"] else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz" if config["use_strelka"] else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz" if config["use_deepvariant"] else None,
        ])),
    conda:
        "../envs/bcftools.yaml"
    params:
        extra="",  # optional
        java_opts="",  # optional
    threads: 56
    resources:
        mem_mb=1024,
    shell:
        """
        echo "##############################################" 
        echo "----   Running BCFTOOLS Concat VCF  ---------" 
        echo "##############################################"

        # Função para processar SNPs e INDELs
        process_vcf() {{
            local snp_file=$1
            local indel_file=$2
            local output_file=$3
            local final_file=$4

            if [ -f "$snp_file" ] || [ -f "$indel_file" ]; then
                echo "Processing $snp_file and $indel_file"
                bcftools concat \
                    --allow-overlaps \
                    "$snp_file" "$indel_file" \
                    -O z \
                    -o "$output_file" \
                    --threads {threads}

                tabix -p vcf "$output_file"

                gatk SelectVariants \
                    --exclude-filtered \
                    -V "$output_file" \
                    -O "$final_file" &
            fi
        }}

        # Contar o número de ferramentas habilitadas
        num_tools=$(echo "{input.snp}" | tr ' ' '\n' | grep -c ".vcf.gz")

        # Processar cada ferramenta habilitada
        for idx in $(seq 0 $((num_tools - 1))); do
            snp_file=$(echo "{input.snp}" | cut -d' ' -f$((idx + 1)))
            indel_file=$(echo "{input.indel}" | cut -d' ' -f$((idx + 1)))
            output_file=$(echo "{output.out}" | cut -d' ' -f$((idx + 1)))
            final_file=$(echo "{output.out_final}" | cut -d' ' -f$((idx + 1)))

            process_vcf "$snp_file" "$indel_file" "$output_file" "$final_file"
        done

        # Esperar todos os processos paralelos terminarem
        wait
        """