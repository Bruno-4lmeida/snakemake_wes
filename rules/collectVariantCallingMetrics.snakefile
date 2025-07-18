rule variantCallingMetrics:
    input:
        vcf=list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz" if config.get("use_deepvariant", False) else None,
        ])),
    output:
        out_vcf=list(filter(None, [
            config["output"] + "output_interomics/metrics/haplotype_caller/{sample}_vcf_metrics.variant_calling_detail_metrics" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/metrics/strelka/{sample}_vcf_metrics.variant_calling_detail_metrics" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/metrics/deep_variant/{sample}_vcf_metrics.variant_calling_detail_metrics" if config.get("use_deepvariant", False) else None,
        ])),
    conda:
        "../envs/picard.yaml"
    threads: 4
    params:
        out_vcf=list(filter(None, [
            config["output"] + "output_interomics/metrics/haplotype_caller/{sample}_vcf_metrics" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/metrics/strelka/{sample}_vcf_metrics" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/metrics/deep_variant/{sample}_vcf_metrics" if config.get("use_deepvariant", False) else None,
        ])),
        ref=config["refgenome"]["ref"],
        intervals=config["refgenome"]["interval"],
        known_var=config["refgenome"]["known_var"],
    shell:
        """
        echo "##############################################" 
        echo "-----  PicardCommandLine CollectVariantCallingMetrics ----" 
        echo "##############################################"

        # Função para processar cada ferramenta
        process_metrics() {{
            local input_vcf=$1
            local output_prefix=$2

            if [ -f "$input_vcf" ]; then
                echo "Processing $input_vcf"
                ### TESTAR AQUI < POIS SO FUNCIONAVA COM O picard.jar > picard CollectVariantCallingMetrics \

                java -jar ~/InterOmics/apps/picard.jar CollectVariantCallingMetrics \
                    I="$input_vcf" \
                    O="$output_prefix" \
                    DBSNP={params.known_var} \
                    TARGET_INTERVALS={params.intervals} &
            fi
        }}

        # Contar o número de ferramentas habilitadas
        num_tools=$(echo "{input.vcf}" | tr ' ' '\n' | grep -c ".vcf.gz")

        # Processar cada ferramenta habilitada
        for idx in $(seq 0 $((num_tools - 1))); do
            input_vcf=$(echo "{input.vcf}" | cut -d' ' -f$((idx + 1)))
            output_prefix=$(echo "{params.out_vcf}" | cut -d' ' -f$((idx + 1)))

            process_metrics "$input_vcf" "$output_prefix"
        done

        # Esperar todos os processos paralelos terminarem
        wait
        """