#!/bin/bash

input_vcf=""
output_vcf=""
sample_name=""

# Lê as flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input_vcf="$2"; shift ;;
        -o|--output) output_vcf="$2"; shift ;;
        -s|--sample) sample_name="$2"; shift ;;
        *) echo "Flag desconhecida: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$input_vcf" || -z "$output_vcf" || -z "$sample_name" ]]; then
    echo "Uso: $0 -i <input_vcf> -o <output_vcf> -s <sample_name>"
    exit 1
fi

# Descompacta se for .vcf.gz
if [[ "$input_vcf" == *.vcf.gz ]]; then
    vcf_stream=$(zcat "$input_vcf")
else
    vcf_stream=$(cat "$input_vcf")
fi

# Cria o novo VCF com cabeçalhos
{
    echo "##fileformat=VCFv4.2"
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name"

    # Processa linhas de variante
    echo "$vcf_stream" | grep -v '^#' | while IFS=$'\t' read -r chrom pos id ref alt qual filter info format sample; do
        # Extrai keys e valores da coluna INFO
        info_keys=$(echo "$info" | tr ';' '\n' | cut -d '=' -f 1 | paste -sd ':' -)
        info_values=$(echo "$info" | tr ';' '\n' | cut -d '=' -f 2 | paste -sd ':' -)

        # Imprime linha formatada
        echo -e "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info_keys\t$info_values\t$sample_name"
    done
} > "$output_vcf"

echo "Arquivo processado: $output_vcf"
