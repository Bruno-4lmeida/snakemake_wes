##exclui env doido
source ~/.bashrc
source ~/.bashrc # elimina envs doidões
#flags do snakemake
# define número de cores
--cores        
#roda o workflow a partir de onde parou refazendo apenas arquivos incompletos 
--rerun-incomplete ou -ri      

# create enviroment
conda env create -f environment.yaml

#controla qual implementação do Conda será usada
--use-conda --conda-frontend conda
--use-conda --conda-frontend mamba

# usa novo config diferente do definido no snakefile
--configfile novo_config.yaml
--config parametro1=valor1 parametro2=valor2
#Não execute nada e exiba o que seria feito.
--dry-reru
# mostra quais regras serão executadas e como seus arquivos de entrada e saída se relacionam (precisa ter instalado o graphviz imagemagick)
snakemake --dag | dot | display

## Util quando adicionar novas amostras ou quando realizar alguma modificacao no workflow
#  a qual nao quer que reanalise os resultados ja existentes
# Execute o comando com --touch uma vez e depois o comando normal

Snakemake --touch # atualiza a data dos arquivos de output
Snakemake --cores all ...

### PARA CONTINUAR A EXECUTAR AS RULES INDEPENDENTES QUANDO OCORRE ERRO COM UMA AMOSTRA

--keep-going




