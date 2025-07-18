rule annovar:
    input:
        list(filter(None, [
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz" if config.get("use_deepvariant", False) else None,
        ])),
    output:
        list(filter(None, [
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annotation_updated.hg38_multianno.csv" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annotation_updated.hg38_multianno.csv" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annotation_updated.hg38_multianno.csv" if config.get("use_deepvariant", False) else None,
        ])),
    threads: config["threads"]["others"]
    params:
        annotation = list(filter(None, [
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annotation_updated" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annotation_updated" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annotation_updated" if config.get("use_deepvariant", False) else None,
        ])),
        out = list(filter(None, [
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annovar.input" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annovar.input" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annovar.input" if config.get("use_deepvariant", False) else None,
        ])),
        info = list(filter(None, [
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annovar.input.INFO" if config.get("use_gatk", False) else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annovar.input.INFO" if config.get("use_strelka", False) else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annovar.input.INFO" if config.get("use_deepvariant", False) else None,
        ])),
        db = "/home/operator/InterOmics/apps/annovar/humandb",
        build_ver = "hg38"
    shell:
        """
        echo "##############################################" 
        echo "-----   ANNOVAR Annotation ------" 
        echo "##############################################"

        # Função para processar cada ferramenta
        process_annovar() {{
            local input_file=$1
            local info_file=$2
            local out_file=$3
            local annotation_file=$4

            if [ -f "$input_file" ]; then
                echo "Processing $input_file"
                /home/operator/InterOmics/apps/annovar/convert2annovar.pl -format vcf4 "$input_file" -includeinfo -outfile "$info_file"
                /home/operator/InterOmics/apps/annovar/convert2annovar.pl -format vcf4 "$input_file" -outfile "$out_file"
                /home/operator/InterOmics/apps/annovar/table_annovar.pl "$out_file" {params.db} -buildver {params.build_ver} -out "$annotation_file" -remove -protocol refGene,cytoBand,wgRna,1000g2015aug_all,gnomad41_genome,exac03,dbscsnv11,dbnsfp47a,clinvar_20240917,cosmic101,avsnp154 -operation g,r,r,f,f,f,f,f,f,f,f -nastring . -csvout &
            fi
        }}

        # Contar o número de ferramentas habilitadas
        num_tools=$(echo "{input}" | tr ' ' '\n' | grep -c ".vcf.gz")

        # Processar cada ferramenta habilitada
        for idx in $(seq 0 $((num_tools - 1))); do
            input_file=$(echo "{input}" | cut -d' ' -f$((idx + 1)))
            info_file=$(echo "{params.info}" | cut -d' ' -f$((idx + 1)))
            out_file=$(echo "{params.out}" | cut -d' ' -f$((idx + 1)))
            annotation_file=$(echo "{params.annotation}" | cut -d' ' -f$((idx + 1)))

            process_annovar "$input_file" "$info_file" "$out_file" "$annotation_file"
        done

        # Esperar todos os processos paralelos terminarem
        wait
        """