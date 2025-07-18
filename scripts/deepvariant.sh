#!/bin/bash

./scripts/deepvariant.sh -l {params.bam_dir} -d {params.ref_dir} -o {params.out_dir} -s {input.bam} -r {input.ref} -t {params.threads} -m {params.model_type} -i {params.bam_index} -g {output.gvcf} -v {output.vcf}

# Função para exibir a ajuda
usage() {
    echo "Uso: $0 -b <bam_dir> -r <ref_dir> -o <output_dir> -s <sample_name> -rn <ref_name> -n <num_shard> -g {output.gvcf} -v {output.vcf}"
    exit 1
}

# Inicializa variáveis
BAM_DIR=""
REF_DIR=""
OUTPUT=""
SAMPLE_NAME=""
REF_NAME=""
NUM_SHARDS=""
BAM_INDEX=""
GVCF=""
VCF=""

# Processa os argumentos
while getopts ":l:d:o:s:r:t:m:i:g:v:" opt; do
    case ${opt} in
        l )
            BAM_DIR=$OPTARG
            ;;
        d )
            REF_DIR=$OPTARG
            ;;
        o )
            OUTPUT=$OPTARG
            ;;
        s )
            SAMPLE_NAME=$OPTARG
            ;;
        r )
            REF_NAME=$OPTARG
            ;;
        t )
            NUM_SHARDS=$OPTARG
            ;;
        m )
            MODEL_TYPE=$OPTARG
            ;;
        i )
            BAM_INDEX=$OPTARG
            ;;    
        g )
            GVCF=$OPTARG
            ;;    
        v )
            VCF=$OPTARG
            ;;    
        \? )
            echo "Opção inválida: -$OPTARG" 1>&2
            echo "opcao invalida"
            usage
            ;;
        : )
            echo "A opção -$OPTARG requer um argumento." 1>&2
            echo "ta faltando um argumento"
            usage
            ;;
    esac
done
shift $((OPTIND -1))

#Verifica se todos os argumentos foram fornecidos
if [ -z "${BAM_DIR}" ] || [ -z "${REF_DIR}" ] || [ -z "${OUTPUT}" ] || [ -z "${SAMPLE_NAME}" ] || [ -z "${REF_NAME}" ] || [ -z "${NUM_SHARDS}" ] || [ -z "${MODEL_TYPE}" ] || [ -z "${BAM_INDEX}" ] || [ -z "${GVCF}" ] || [ -z "${VCF}" ]; then
    echo "Nem todos argumentos foram fornecidos"
    usage
fi

mkdir -p /home/operator/InterOmics/workflows/interomics_snakemake/output_interomics/vcf/call/
# Usa as variáveis
echo "Diretorio arquivo BAM: $BAM_DIR"
echo "Diretorio Arquivo de Referência: $REF_DIR"
echo "Diretório de Saída (output): $OUTPUT"
echo "Nome da Amostra (path + sample): $SAMPLE_NAME"
echo "Nome do genoma de referencia: $REF_NAME"
echo "Número de Shards: $NUM_SHARDS"
echo "Model type: $MODEL_TYPE"
echo "index do arquivo BAM: $BAM_INDEX"
echo "Output gvcf: $GVCF"
echo "Output vcf: $VCF"


# RUN:
docker run \
     --rm \
     -v ${BAM_DIR}:${BAM_DIR} \
     -v ${REF_DIR}:${REF_DIR} \
     -v ${OUTPUT}:output/vcf/ \
     -v ${BAM_INDEX}:${BAM_INDEX} \
     google/deepvariant:latest \
     /opt/deepvariant/bin/run_deepvariant \
     --num_shards=${NUM_SHARDS} \
     --model_type=${MODEL_TYPE} \
     --ref=${REF_NAME} \
     --reads=${SAMPLE_NAME} \
     --output_vcf=output/vcf/COL53_subset_reads.vcf.gz \
     --output_gvcf=output/vcf/COL53_subset_reads.g.vcf.gz
