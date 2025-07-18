library(tidyverse)
library(data.table)
library(DBI)
library(RPostgres)


#TSV file
for(i in 1:length(
  list.files(pattern = "hg38_multianno.csv"))) {
  assign(paste0("file",i),
         fread(list.files(pattern = "hg38_multianno.csv")[i]) %>% 
           mutate(sample = str_split_1(list.files(pattern = "hg38_multianno.csv")[i], "_")[1])
  )
  
}
rm(i)

dataset_csv = eval(parse(text = 
                           paste0("bind_rows(",paste0(ls(pattern = "^file"), collapse = ","),")")
))
rm(list = ls(pattern = "^file"))

#NORMAL avinput file
file_list = list.files(
  pattern = "annovar.input")[!grepl("INFO|csv", list.files(pattern = "annovar.input"))]

for(i in 1:length(file_list)) {
  
  
  assign(paste0("file",i),
         fread(file_list[i]) %>% 
           mutate(sample = str_split_1(list.files(pattern = "annovar.input")[i], "_")[1]) %>% 
           # mutate(snv_indel = str_split_1(list.files(pattern = "annovar.input")[i], "_")[3]) %>% 
           # rowwise() %>% 
           # mutate(snv_indel = str_split_1(snv_indel, "\\.")[1]) %>% 
           dplyr::select(genotype = V6)
  )         
}
rm(i, x, file_list)

dataset_input = eval(parse(text = 
                             paste0("bind_rows(",paste0(ls(pattern = "^file"), collapse = ","),")")
))
rm(list = ls(pattern = "^file"))

#INFO NORMAL avinput file
file_list = list.files(
  pattern = "_annovar.input.INFO")[grepl("INFO", list.files(pattern = "_annovar.input.INFO"))]

for(i in 1:length(file_list)) {
  
  assign(paste0("file",i),
         fread(file_list[i]) %>% 
           mutate(sample = str_split_1(list.files(pattern = "_annovar.input.INFO")[i], "_")[1]) %>%
           # mutate(snv_indel = str_split_1(list.files(pattern = "_annovar.input.INFO")[i], "_")[3]) %>%
           # rowwise() %>%
           # mutate(snv_indel = str_split_1(snv_indel, "\\.")[1]) %>% 
           dplyr::select(V12:V15) %>% 
           dplyr::rename(filter = V12) %>%
           dplyr::rename(format = V13) %>% 
           dplyr::rename(info_header = V14) %>% 
           dplyr::rename(info = V15)
         
  )
}

rm(i, x, file_list)

dataset_info_input = eval(parse(text = 
                                  paste0("bind_rows(",paste0(ls(pattern = "^file"), collapse = ","),")")
))
rm(list = ls(pattern = "^file"))

####clean dataset####
dataset = bind_cols(dataset_csv, 
                    dataset_input, 
                    dataset_info_input) %>% 
  janitor::clean_names() %>% 
  mutate(across(where(is.character), 
                ~ na_if(.x, "."))) %>% 
  mutate(across(gnomad41_genome_af:gnomad41_genome_af_sas, as.numeric)) %>% 
  mutate(across(ends_with("score"), as.numeric)) %>% 
  mutate(cadd_phred = as.numeric(cadd_phred)) %>% 
  dplyr::rename(start_bp = start) %>% 
  dplyr::rename(end_bp = end) %>% 
  mutate(avsnp154 = coalesce(avsnp154, 
                             str_c(chr, start_bp, end_bp, ref, alt, sep = ":")))

rm(dataset_csv, dataset_input, dataset_info_input)

#function to organize INFO/FORMAT

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

#separate INFO/FORMAT colums

format_to_separate = dataset %>% 
  rowwise() %>% 
  mutate(to_separate = reconstruct_format(info_header, info)) %>% 
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


dataset = dataset %>% 
  rownames_to_column("index") %>%
  left_join(., format_to_separate, join_by(index == index)) %>%
  select(-c(index, info_header, info)) %>%
  mutate_at(.vars = c("DP", "GQ", "GQX", "SB"),
            .funs = ~ as.numeric(.x))

rm(format_to_separate, reconstruct_format)

# write table

write.table(dataset,
            file = "dataset.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

#### load clean dataframe to database ####

#create a connection to postgreSQL database
con = dbConnect(drv = RPostgres::Postgres(),
                dbname="liquid_biopsy",
                host = "10.23.5.101",
                port = 5432,
                user = "nour",
                password = "nour")

# Write the R dataframe to the PostgreSQL table
dbWriteTable(conn = con,
             name = "somatic_calling",   
             value = dataset,
             overwrite = TRUE)

# Close the connection
dbDisconnect(con)
rm(con)


