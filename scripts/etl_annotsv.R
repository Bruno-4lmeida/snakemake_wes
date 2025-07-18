#load libraries and functions
library(tidyverse)

reconstruct_format = function(info_header, info){
  
  # Load the stringr package
  require(stringr)
  
  # Define the two strings
  keys = info_header
  values = info
  
  # Split the strings by ":"
  keys_split = str_split(keys, ":", simplify = TRUE)
  values_split = str_split(values, ":", simplify = TRUE)
  
  # Merge the corresponding elements
  return(paste0(keys_split, "=", values_split, collapse = ";"))
  
}

#load annotation files
for(i in list.files(pattern = "_annotsv_SVs.tsv")){
  assign(str_replace(i, "_annotsv_SVs.tsv", ""),
         read_delim(i) %>% 
           mutate(sample = str_replace(i, "_annotsv_SVs.tsv", "")) %>% 
           janitor::clean_names() %>% 
           mutate(sv_chrom = paste0("chr", sv_chrom)) %>% 
           rename(format_header = format,
                  format = 15) %>% 
           select(annot_sv_id, sample, id, 
                  sv_chrom, sv_start, sv_end, sv_length,
                  cyto_band, gene_name, gene_count,
                  frameshift, location, location2,
                  gnom_ad_p_li,
                  gen_cc_disease, gen_cc_classification, gen_cc_pmid,
                  omim_id, omim_phenotype, omim_inheritance, 
                  acmg, acmg_class,
                  starts_with("p_"),
                  info, format_header, format)
         )
}
rm(i)

dataset = eval(parse(text = paste0(
  "bind_rows(", paste0(setdiff(ls(), "reconstruct_format"), collapse = ","), ")")))


Cholestasis=bind_rows(COL12,COL14,COL16,COL18,COL20,COL21,COL22,COL23,
                    COL24,COL3,COL3637M,COL3637P,COL39M,COL39P,col52,
                    COL8,COL9,mwes769,mwes770,mwes771,mwes772,mwes773
                    ,wes83,wes84,wescol47,wescol47M,wescol47P)
Cholestasis=Cholestasis %>%
  filter(str_detect(format, "PASS")) %>%
  filter(str_detect(gen_cc_disease, "cholestasis"))

### para o manta
padroniza_omim <- function(df) {
  df$omim_id <- as.character(df$omim_id)
  return(df)
}


Cholestasis_manta <- bind_rows(
  padroniza_omim(COL12), padroniza_omim(COL14), padroniza_omim(COL16),
  padroniza_omim(COL18), padroniza_omim(COL20), padroniza_omim(COL21),
  padroniza_omim(COL22), padroniza_omim(COL23), padroniza_omim(COL24),
  padroniza_omim(COL3), padroniza_omim(COL3637M), padroniza_omim(COL3637P),
  padroniza_omim(COL39M), padroniza_omim(COL39P), padroniza_omim(col52),
  padroniza_omim(COL8), padroniza_omim(COL9), padroniza_omim(mwes769),
  padroniza_omim(mwes770), padroniza_omim(mwes771), padroniza_omim(mwes772),
  padroniza_omim(mwes773), padroniza_omim(wes83), padroniza_omim(wes84),
  padroniza_omim(wescol47), padroniza_omim(wescol47M), padroniza_omim(wescol47P)
)


Cholestasis_manta=Cholestasis_manta %>%
  filter(str_detect(gen_cc_disease, "cholestasis"))  ## DO NOT USE PASS FILTER ON MANTA

### MANTA AND DELLY COMBINED
Cholestasis_combined <- bind_rows(Cholestasis, Cholestasis_manta) %>%
  distinct(sample, gene_name, .keep_all = TRUE)

library(dplyr)

Cholestasis_combined <- semi_join(
  Cholestasis,
  Cholestasis_manta,
  by = c("sample", "gene_name")
)





#ETL
dataset_filtered = dataset %>% 
  separate(info, 
           into = c("IMPRECISE", "SVTYPE", "SVMETHOD", "END", "CIPOS", "CIEND", "MP"),
           sep = ";") %>% 
  mutate_at(.vars = c("IMPRECISE", "SVTYPE", "SVMETHOD", "END", "CIPOS", "CIEND", "MP"),
            .funs = ~ str_replace_all(.x, 
                                     "IMPRECISE|SVTYPE|SVMETHOD|END|CIPOS|CIEND|MP|=", "")) %>% 
  janitor::clean_names()

format_to_separate = dataset_filtered %>% 
  rowwise() %>% 
  mutate(to_separate = reconstruct_format(format_header, format)) %>% 
  select(to_separate) %>% 
  rownames_to_column("index") %>% 
  mutate(to_separate = str_squish(to_separate)) %>% 
  mutate(to_separate = str_replace(to_separate, ";$", "")) %>% 
  mutate(to_separate = str_replace_all(to_separate, "\"", "")) %>% 
  mutate(to_separate = str_replace_all(to_separate, "; ", ";")) %>% 
  separate_rows(to_separate, sep = ";") %>%
  separate(to_separate, 
           into = c("key", "value"), 
           sep = "=", 
           fill = "right") %>%
  pivot_wider(names_from = key,
              values_from = value, 
              values_fill = NA,
              values_fn = 
                list(value = function(x) paste(unique(x), collapse = ",")))

dataset_filtered = dataset_filtered %>% 
  rownames_to_column("index") %>%
  left_join(., format_to_separate, join_by(index == index)) %>%
  select(-c(index, format_header, format)) %>% 
  janitor::clean_names()

dataset_pass = dataset_filtered %>%
  filter(ft == "PASS") %>% 
  mutate(cyto_band = paste0(sv_chrom, cyto_band))


dataset_pass %>%
  write.table(., file = "dataset_SVs_pass.tsv", sep = "\t", 
              row.names = FALSE, quote = FALSE, col.names = TRUE)

