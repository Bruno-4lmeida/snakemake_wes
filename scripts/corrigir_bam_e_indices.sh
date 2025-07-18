#!/bin/bash

# Diretório com os arquivos BAM
BAM_DIR="."

# Função para verificar o cabeçalho do arquivo BAM
check_chr_prefix() {
    local bam_file="$1"
    samtools view -H "$bam_file" | grep -q "SN:chr"
}

# Processar cada arquivo BAM no diretório
for bam_file in "$BAM_DIR"/*.bam; do
    if [[ -f "$bam_file" ]]; then
        echo "Verificando o arquivo: $bam_file"

        # Extrair cabeçalho para diagnóstico
        header_content=$(samtools view -H "$bam_file")
        
        # Registrar cabeçalho em um arquivo para análise posterior
        echo "$header_content" > "${bam_file%.bam}_header.txt"

        if check_chr_prefix "$bam_file"; then
            echo "O arquivo $bam_file contém o prefixo 'chr'. Corrigindo..."

            # Extrai o cabeçalho, remove o prefixo "chr", e reinsere o cabeçalho corrigido
            echo "$header_content" | sed 's/SN:chr\([0-9XYMT]\)/SN:\1/' > header.sam
            samtools reheader header.sam "$bam_file" > "${bam_file%.bam}_corrigido.bam"
            
            # Substitui o arquivo original pelo corrigido
            mv "${bam_file%.bam}_corrigido.bam" "$bam_file"
            
            # Limpa o cabeçalho temporário
            rm header.sam
            
            # Recria os arquivos de índice (.bai e .sbi)
            echo "Recriando os índices para $bam_file..."
            samtools index "$bam_file"  # Gera o arquivo .bai
            samtools index -c "$bam_file"  # Gera o arquivo .sbi (caso necessário)
            
            echo "Correção e reindexação concluídas para: $bam_file"
        else
            echo "O arquivo $bam_file já está no formato numérico. Nenhuma ação necessária."
        fi
    fi
done

echo "Processamento concluído."
