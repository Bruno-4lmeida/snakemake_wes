# snakemake_wes — Workflow Snakemake para análise de Exoma Completo (WES)

Este repositório contém um workflow em **Snakemake** para processar dados de **Whole Exome Sequencing (WES)** a partir de arquivos **FASTQ pareados**, produzir BAMs processados (alinhados, com duplicatas marcadas e recalibrados) e executar **variant calling** (germline ou somático) com etapas de **seleção, filtragem, concatenação e anotação**, além de um conjunto de **métricas de qualidade** e um relatório agregador com **MultiQC**.

O pipeline é configurável via `config.yaml`, permitindo alternar **alinhador** (BWA ou DRAGMAP) e **chamador de variantes** (GATK HaplotypeCaller, Strelka, ou DeepVariant), além de escolher entre análise **germline** ou **somática**.

---

## Visão geral do workflow

O `Snakefile` principal:

- Lê as configurações do `config.yaml`
- Detecta automaticamente as amostras a partir do diretório de entrada
- Inclui as rules necessárias de acordo com as opções selecionadas

Fluxo geral:

1. **Detecção de amostras**
2. **Pré-processamento e QC** (fastp + MultiQC)
3. **Alinhamento** (BWA ou DRAGMAP)
4. **Pós-processamento do BAM** (MarkDuplicates + BQSR)
5. **Variant calling** (germline: GATK/Strelka/DeepVariant; somático: Strelka somático)
6. **Pós-processamento de VCF** (seleção, filtragem, concatenação)
7. **Anotação** (ANNOVAR)
8. **Métricas** (GC/HS/WGS + métricas de variant calling) e consolidação em CSV

---

## Estrutura do repositório (módulos principais)

O workflow é modular e inclui diferentes arquivos em `rules/` dependendo do `config.yaml`.

### Seleção do alinhador

- `use_bwa: True`  → `rules/bwa_mem.snakefile`
- `use_bwa: False` → `rules/dragmap.snakefile`

### Germline vs Somático

- `germline: True`:
  - `use_gatk: True` → `rules/haplotypeCaller.snakefile`
  - `use_strelka: True` → `rules/strelka.snakefile`
  - `use_deepvariant: True` → `rules/deepvariant_singularity.snakefile` (somente se Singularity estiver instalado)
- `germline: False`:
  - `rules/strelka_somatic.snakefile`

### Rules sempre incluídas

Independente do caller, o workflow sempre inclui:

- `rules/multiqc.snakefile`
- `rules/fastp.snakefile`
- `rules/markduplicates.snakefile`
- `rules/recal_table.snakefile`
- `rules/recal_apply.snakefile`
- `rules/selectVar.snakefile`
- `rules/filterVar.snakefile`
- `rules/concat_bcftools.snakefile`
- `rules/annovar.snakefile`
- `rules/collectGCBiasMetrics.snakefile`
- `rules/collectHSmetrics.snakefile`
- `rules/collectWGSMetrics.snakefile`
- `rules/collectVariantCallingMetrics.snakefile`
- `rules/process_wgs_metricsR.snakefile`
- `rules/process_gc_metrics.snakefile`
- `rules/process_hs_metrics.snakefile`

---

## Requisitos

### Dependências gerais

- **Snakemake**
- **Conda/Mamba** (recomendado) caso você use ambientes (por regra ou global)
- Ferramentas utilizadas pelas rules (ex.: `fastp`, `bwa`/`dragmap`, `gatk`, `bcftools`, `annovar`, ferramentas de métricas etc.)

> Observação: as dependências exatas dependem do conteúdo das rules e/ou de ambientes Conda definidos no repositório.

### DeepVariant (opcional)

Se `use_deepvariant: True`, o workflow verifica se o **Singularity** está instalado e acessível via:

- `singularity --version`

Se estiver instalado, usa o `.sif` definido em `deepvariant_sif_path`.  
Caso contrário, a execução do DeepVariant não é possível e você deve instalar Singularity ou escolher outro chamador.

---

## Entradas (input) e detecção automática de amostras

O pipeline infere amostras automaticamente usando:

- Diretório de entrada: `config["input"]`
- Padrão obrigatório para R1: `{sample}_1.fastq.gz`
- O pipeline assume que existe também o pareado: `{sample}_2.fastq.gz`

Exemplo de nomes válidos no diretório de entrada:

- `SAMPLE01_1.fastq.gz`
- `SAMPLE01_2.fastq.gz`

Qualquer arquivo que case com `*_1.fastq.gz` define uma amostra e entra no DAG. O `Snakefile` imprime a lista de amostras detectadas ao iniciar.

---

## Configuração (`config.yaml`)

A execução do pipeline é controlada por `config.yaml`.

### Modo de análise

- `germline: True` → fluxo germline (com GATK/Strelka/DeepVariant conforme flags)
- `germline: False` → fluxo somático (Strelka somático)

### Alinhamento

- `use_bwa: True` → usa BWA
- `use_bwa: False` → usa DRAGMAP

### Variant calling (germline)

- `use_gatk: True` → GATK HaplotypeCaller
- `use_strelka: True` → Strelka (germline)
- `use_deepvariant: True` → DeepVariant (via Singularity)

> Você pode habilitar mais de um caller ao mesmo tempo. O `rule all` exige os outputs de todos os callers habilitados.

### DeepVariant: tipo de modelo

- `model_type: WES` (usado quando DeepVariant está habilitado)

### Caminhos gerais

> **Importante:** os paths foram definidos para terminar com `/` no YAML.

- `input: /home/operator/InterOmics/workflows/interomics_snakemake/samples/`
- `output: /home/operator/InterOmics/workflows/interomics_snakemake/`

### Referência e recursos

Bloco `refgenome`:

- `ref_dir`: diretório de referência
- `ref`: FASTA de referência (GRCh38)
- `bed`: BED do kit de captura do exoma
- `interval`: interval list (GATK)
- `dict`: dicionário do FASTA
- `known_var`: VCF de variantes conhecidas (dbSNP) usado na BQSR
- `data_source`: recursos adicionais (ex.: Funcotator germline), conforme rules

### Recursos de execução (threads e memória)

Exemplo do YAML fornecido:

- `threads`:
  - `others: 1`
  - `bwa: 13`
  - `deepvariant: 56`
  - `recal_table: 8`
  - `recal_apply: 8`
  - (`gatk` aparece no YAML, mas o valor não foi incluído no trecho enviado)

- `ram_memory`:
  - `markduplicates: 102400`
  - `recal_table: 102400`
  - (`bwa` e `gatk` aparecem no YAML, mas os valores não foram incluídos no trecho enviado)

> Recomenda-se completar os campos ausentes (`threads.gatk`, `ram_memory.bwa`, `ram_memory.gatk`) caso suas rules façam uso desses parâmetros.

---

## Como executar

Execute a partir da raiz do repositório (onde estão `Snakefile` e `config.yaml`).

### Execução local (exemplo)

```bash
snakemake --cores 16 --use-conda --rerun-incomplete --printshellcmds
