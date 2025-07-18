docker run \
     --rm \
     -v  "/home/operator/InterOmics/workflows/interomics_snakemake/output_interomics/recal_apply":"/input" \
     -v "/home/operator/InterOmics/workflows/interomics_snakemake/reference":"/reference" \
     -v "/home/operator/InterOmics/workflows/interomics_snakemake/output_interomics/vcf/call":"/output" \
     google/deepvariant:latest \
     /opt/deepvariant/bin/run_deepvariant \
     --num_shards=50 \
     --model_type=WES \
     --ref=/reference/GRCh38_numeric_genome.fasta \
     --reads=/input/COL53_subset_reads_markduplicates_recal.bam \
     --output_vcf=/output/COL53_subset_reads.vcf.gz \
     --output_gvcf=/output/COL53_subset_reads.g.vcf.gz

     #####################################

for i in $(ls | grep _markduplicates_recal.bam | sed 's/_markduplicates_recal.bam//' | uniq);
do


ulimit -u 10000 
BIN_VERSION="latest"

SAMPLE=$i

docker run \
-v "/home/operator/InterOmics/workflows/interomics_snakemake/reference/":"/db" \
-v "/home/operator/Downloads/ab_analitica/output_interomics/recal_apply/":"/input" \
-v "/home/operator/Downloads/ab_analitica/output_interomics/vcf/call:/output" \
google/deepvariant:"${BIN_VERSION}" \
/opt/deepvariant/bin/run_deepvariant \
--num_shards=50 \
--model_type=WES \
--ref=/db/GRCh38_numeric_genome.fasta \
--reads=/input/"${SAMPLE}"_markduplicates_recal.bam \
--regions=/db/hg38_no_alt.bed \
--output_vcf=/output/"${SAMPLE}".vcf.gz 


docker stop $(docker ps -a -q)
docker rm $(docker ps -a -q)

####make examples
# real    11m55.133s
# user    554m45.872s
# sys     4m12.451s

####call variants
# real    19m11.071s
# user    736m7.289s
# sys     4m3.482s

####post processing
# real    1m8.630s
# user    2m51.279s
# sys     0m28.120s

done


