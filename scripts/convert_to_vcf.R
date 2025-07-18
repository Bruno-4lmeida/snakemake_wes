# Carregar bibliotecas necessárias
library(dplyr)
library(readr)

# Caminhos dos arquivos
input_csv <- "/home/operator/InterOmics/workflows/interomics_snakemake/output_interomics/cnvs/cnmops/resultscnmops_cnvs.csv"
output_vcf <- "/home/operator/InterOmics/workflows/interomics_snakemake/output_interomics/cnvs/cnmops/resultscnmops_cnvs.vcf"

# Ler o CSV
cnvs <- read_csv(input_csv)

# Criar cabeçalho VCF
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">",
  "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">",
  "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

# Converter para formato VCF
vcf_body <- cnvs %>%
  mutate(
    CHROM = seqnames,
    POS = start,
    ID = ".",
    REF = "N",
    ALT = ifelse(median > 0, "<DUP>", "<DEL>"),
    QUAL = ".",
    FILTER = "PASS",
    INFO = paste0(
      "SVLEN=", end - start, ";",
      "SVTYPE=", ifelse(median > 0, "DUP", "DEL"), ";",
      "CN=", median, ";",
      "SAMPLE=", sampleName
    )
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
  arrange(as.numeric(CHROM), POS)  # Ordenar por cromossomo e posição

# Escrever VCF
writeLines(vcf_header, output_vcf)
write.table(vcf_body, output_vcf, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

cat("Arquivo VCF gerado com sucesso:", output_vcf, "\n")
