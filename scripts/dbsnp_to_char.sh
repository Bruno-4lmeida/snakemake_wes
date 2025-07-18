#!/bin/bash

# Caminho para o arquivo original e o novo arquivo
input_file="reference/dbsnp_156.vcf.gz"
output_file="reference/dbsnp_156_chr.vcf.gz"

# Descompacta, ajusta os cromossomos e recompacata
zcat "$input_file" | \
sed -E 's/(##contig=<ID=)([0-9]+)(>)/\1chr\2\3/; s/^([0-9]+)\b/chr&/' | \
bgzip > "$output_file"

# Indexa o novo arquivo VCF (opcional)
tabix -p vcf "$output_file"
