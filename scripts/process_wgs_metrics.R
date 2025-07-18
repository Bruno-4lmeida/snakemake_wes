library(tidyverse)

# Lê o argumento da linha de comando para o diretório base
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Erro: Nenhum diretório de entrada foi fornecido como argumento.")
}
base_dir = args[1]
base_dir = gsub("/$", "", base_dir)  # Garante que o diretório termine com uma barra

out_dir = args[2]

# Inicializa um data.frame vazio
ngs_metrics = data.frame()

# Lista os arquivos de métricas
file_list = list.files(path = base_dir, pattern = "_collect_wgs_metrics.txt$")

# Valida se há arquivos
if (length(file_list) == 0) {
  stop("Nenhum arquivo correspondente foi encontrado no diretório especificado.")
}

# Processa os arquivos
for (i in file_list) {
  x = read_delim(
    paste0(base_dir,"/", i),
    skip = 6,
    n_max = 1,
    delim = "\t",
    locale = locale(decimal_mark = ","),
    show_col_types = FALSE
  )
  print(paste0(base_dir, i))  # Para verificar o caminho

  
  # Adiciona o nome do arquivo
  x = x %>%
    mutate(filename = i) %>%
    relocate(filename, .before = everything())
  
  # Combina os resultados
  ngs_metrics = bind_rows(ngs_metrics, x)
}

# Caminho de saída
output_file = paste0(out_dir)

# Salva o resultado
write_csv(ngs_metrics, output_file)
