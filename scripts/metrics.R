library(tidyverse)

base_dir = "~/Downloads/nextseq/output_interomics/metrics/"

#### NGS metrics ####
ngs_metrics = NULL

for(i in list.files(path = base_dir,
                    pattern = "_collect_wgs_metrics.txt$")){
  x = read_delim(paste0(base_dir, i), 
                 skip = 6,
                 n_max = 1,
                 delim = "\t", 
                 # locale = locale(decimal_mark = ","),
                 show_col_types = FALSE)
  x = x %>% 
    mutate(filename = i) %>% 
    relocate(filename, .before = GENOME_TERRITORY)
  
  ngs_metrics = bind_rows(ngs_metrics, x)
  
  rm(x)
}

ngs_metrics %>% 
  janitor::clean_names() %>% 
  write.table(.,
              file = paste0(base_dir, "/summary/collect_wgs_metrics.csv"), 
              quote = FALSE, 
              row.names = FALSE, 
              sep = ";",
              dec = ".")

#### GC metrics ####
gc_metrics = NULL

for(i in list.files(path = base_dir,
                    pattern = "_gc_summary_metrics.txt$")){
  x = read_delim(paste0(base_dir, i), 
                 skip = 6,
                 n_max = 1,
                 delim = "\t", 
                 locale = locale(decimal_mark = ","),
                 show_col_types = FALSE)
  x = x %>% 
    mutate(filename = i) %>% 
    relocate(filename, .before = ACCUMULATION_LEVEL)
  
  gc_metrics = bind_rows(gc_metrics, x)
  
  rm(x)
}

gc_metrics %>% 
  janitor::clean_names() %>% 
  write.table(.,
              file = paste0(base_dir, "/summary/gc_metrics.csv"), 
              quote = FALSE, 
              row.names = FALSE, 
              sep = ";",
              dec = ".")

#### HS metrics ####

hs_metrics = NULL

for(i in list.files(path = base_dir,
                    pattern = "_hs_metrics.txt$")){
  x = read_delim(paste0(base_dir, i), 
                 skip = 6,
                 n_max = 1,
                 delim = "\t", 
                 # locale = locale(decimal_mark = ","),
                 show_col_types = FALSE)
  
  x = x %>% 
    mutate(filename = i) %>% 
    relocate(filename, .before = 1) %>% 
    mutate(PCT_TARGET_BASES_1000X = as.character(PCT_TARGET_BASES_1000X))
  
  hs_metrics = bind_rows(hs_metrics, x)
  
  rm(x)
}


hs_metrics %>% 
  janitor::clean_names() %>% 
  write.table(.,
              file = paste0(base_dir, "/summary/hs_metrics.csv"), 
              quote = FALSE, 
              row.names = FALSE, 
              sep = ";",
              dec = ".")

####proportions####

prop_ad10 = read_delim(paste0(base_dir, "proportions.csv"))

prop_ad10 = prop_ad10 %>% 
  mutate(Proportion = str_replace(Proportion, 
                                  "Proportion of variants with FILTER=PASS and DP>=10: ",
                                  "0")) %>% 
  mutate(Proportion = as.numeric(Proportion))

prop_ad10 %>% 
  janitor::clean_names() %>% 
  write.table(.,
              file = paste0(base_dir, "/summary/prop_ad10_bed.csv"), 
              quote = FALSE, 
              row.names = FALSE, 
              sep = ";",
              dec = ".")

####proportions bed ####

prop_ad10_bed = read_delim("~/Downloads/nextseq/output_interomics/metrics/proportions.csv")

prop_ad10_bed = prop_ad10_bed %>% 
  mutate(Proportion = str_replace(Proportion, 
                                  "Proportion of variants with FILTER=PASS and DP>=10: ",
                                  "0")) %>% 
  mutate(Proportion = as.numeric(Proportion))

prop_ad10_bed %>% 
  janitor::clean_names() %>% 
  write.table(.,
              file = paste0(base_dir, "/summary/prop_ad10_bed.csv"), 
              quote = FALSE, 
              row.names = FALSE, 
              sep = ";",
              dec = ".")


#### coverage uniformity ####

# Função para calcular a uniformidade de cobertura para uma única amostra
calculate_uniformity <- function(sample_name, coverage_dir, interval_file) {
  
  # Definir o caminho do arquivo de cobertura da amostra
  coverage_file <- file.path(coverage_dir, paste0(sample_name, "_per_base_coverage.txt.gz"))
  
  # Verificar se os arquivos existem
  if (!file.exists(coverage_file) || !file.exists(interval_file)) {
    stop(paste("Erro: Arquivo(s) não encontrado(s) para a amostra", sample_name))
  }
  
  # Ler arquivos e calcular métricas com tidyverse
  coverage_data <- fread(coverage_file)
  
  intervals <- read_delim(interval_file, skip = 196, col_names = FALSE)
  
  # Filtrar apenas as bases-alvo no exoma e calcular métricas
  results <- coverage_data %>%
    inner_join(intervals, join_by(chrom == X1, pos >= X2, pos <= X3)) %>%
    summarise(
      Total_Bases = n(),
      Mean_Coverage = mean(coverage),
      Threshold = Mean_Coverage * 0.2,
      Bases_Above_Threshold = sum(coverage > Threshold),
      Uniformity_Coverage = (Bases_Above_Threshold / Total_Bases) * 100
    ) %>%
    mutate(Sample = sample_name) %>%
    dplyr::select(Sample, Threshold, Total_Bases, Bases_Above_Threshold, Uniformity_Coverage) %>%
    janitor::clean_names()
  
  return(results)
}

coverage_dir <- "~/Downloads/nextseq/output_interomics/metrics/"  # Diretório onde estão os arquivos das amostras
interval_file <- "/home/operator/InterOmics/workflows/interomics_snakemake/reference/hg38_Twist_ILMN_Exome_2.5_Panel.interval_list"  # Caminho para o arquivo de intervalos
output_file <- "uniformity_coverage_results.csv"

sample_names = list.files(coverage_dir, 
                          pattern = "_per_base_coverage.txt.gz$") %>%
  str_remove("_per_base_coverage.txt.gz")

uniformity = NULL
for(i in sample_names){
  x = calculate_uniformity(sample_name = i, coverage_dir = coverage_dir, interval_file = interval_file)
  uniformity = bind_rows(x, uniformity)
  rm(x)
}

uniformity %>% 
  write.table(., 
              file = paste0(base_dir, "summary/coverage_uniformity_bed.csv"),
              sep = ";", 
              
              quote = FALSE, 
              row.names = FALSE,
              col.names = TRUE)
