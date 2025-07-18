# run_cnmops.R

# Carregar pacotes necessários
library(cn.mops)
library(Rsamtools)
library(parallel)

# Receber argumentos da linha de comando
args <- commandArgs(trailingOnly = TRUE)

bam_dir <- as.character(args[1])
output_file <- args[2]
threads <- args[3]
bed_file <- args[4]

setwd(bam_dir) # set dir

if (!file.exists(bam_dir)) stop("The given directory was not founding.")

# BAMFIles
BAMFiles <- list.files(path =bam_dir ,pattern=".bam$")

# Print on terminal
cat("BAM files: ", BAMFiles, "\n")

# read bed file
segments <-  read.csv2(bed_file, header=FALSE, sep="")
segments <- GRanges(segments[,1],IRanges(as.numeric(segments[,2]),as.numeric(segments[,3])))

# reading the segments from bam files using the coordinates in bed
input_data <- getSegmentReadCountsFromBAM(BAMFiles,GR=segments, parallel = TRUE)

# running exome cn.mops
cnv_result <- exomecn.mops(input_data, parallel = as.numeric(threads))
cnv_result <- calcIntegerCopyNumbers(cnv_result)
cnvs <- cnv_result %>%
  as.data.frame() %>% 
  mutate(CN = gsub("CN", "", CN)) %>% 
  arrange(as.numeric(gsub("chr", "", seqnames))) 

cnvr=as.data.frame(cnvr(cnv_result))

### FOR WHOLE GENOME SEQUENCING ####
# Select chromossomes
#chromossomes <- c("chr1", "chr2", "chr3") #Mudar aqui tambem
#chromossomos <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
#                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
#                  "chr20", "chr21", "chr22")
#input_data=getReadCountsFromBAM(BAMFiles[1:6], refSeqNames=chromossomes, parallel = TRUE) ##mudar aqui remover de 1:6
# Configurar para usar a versão paralela
#cnv_result <- cn.mops(input = input_data, parallel = as.numeric(threads))
#cnv_result <- calcIntegerCopyNumbers(cnv_result)

#### FOR HAPLOID GENOMES ChrX and ChrY
#resHaplo <- haplocn.mops(X)
#resHaplo <- calcIntegerCopyNumbers(resHaplo)




##### CONVERT TO VCF FORMAT
# Carregar bibliotecas
library(dplyr)

library(dplyr)

convert_to_vcf <- function(cnmops_results, output_file) {
  # Cabeçalho do VCF
  vcf_header <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">",
    "##INFO=<ID=CN,Number=1,Type=Integer,Description=\"Copy number\">",
    "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  )
  
  # Criar o corpo do VCF
  vcf_body <- cnmops_results %>%
    mutate(
      CHROM = seqnames,  
      POS = start,       
      ID = ".",          
      REF = "N",         
      ALT = ifelse(CN > 2, "<DUP>", "<DEL>"),  
      QUAL = ".",        
      FILTER = "PASS",   
      INFO = paste0(
        "SVLEN=", end - start, ";",
        "SVTYPE=", ifelse(CN > 2, "DUP", "DEL"), ";",
        "CN=", CN, ";",
        "SAMPLE=", sampleName  # Adicionando o nome da amostra ao INFO
      )
    ) %>%
    select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
  
  # Salvar o arquivo VCF
  writeLines(vcf_header, output_file)
  write.table(vcf_body, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

# OUTPUT VCF FILES
convert_to_vcf(cnvs, gsub(x=output_file, "_cnvs.csv", "_cnvs.vcf"))

# OUTPUT: global CNVs e CNVRs
write.csv(cnvs, file = output_file, row.names = FALSE)
write.csv(cnvr, file = gsub(x=output_file, "_cnvs", "_cnvr"), row.names = FALSE)
