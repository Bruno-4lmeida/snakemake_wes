export home_dir=~/InterOmics
export picard=${home_dir}/apps
export db=${home_dir}/workflows/interomics_snakemake/reference

export metrics_dir=../QC
export scripts=~/InterOmics/workflows/interomics_snakemake/scripts

java -Xmx10G -jar ${picard}/picard.jar BedToIntervalList \
      I=/home/operator/InterOmics/workflows/interomics_snakemake/reference/hg38_no_alt.bed.gz \
      O=/home/operator/InterOmics/workflows/interomics_snakemake/reference/hg38_no_alt.interval_list \
      SD=${db}/GRCh38_character_genome.dict

#NGS metrics
ls | grep _recal.bam | sed 's/_recal.bam//' | parallel --max-args=1 \
java -Xmx10G -jar ${picard}/picard.jar CollectWgsMetrics \
      I={1}_recal.bam \
      O=${metrics_dir}/{1}_collect_wgs_metrics.txt \
      R=${db}/GRCh38_character_genome.fasta \
      INTERVALS=${db}/hg38_no_alt.interval_list \
      MINIMUM_MAPPING_QUALITY=20 \
      MINIMUM_BASE_QUALITY=20 \
      COVERAGE_CAP=250

#GC metrics
ls | grep _recal.bam | sed 's/_recal.bam//' | parallel --max-args=1 \
java -Xmx8G -jar ${picard}/picard.jar CollectGcBiasMetrics \
      I={1}_recal.bam \
      R=${db}/GRCh38_character_genome.fasta \
      O=${metrics_dir}/{1}_gc_bias_metrics.txt \
      CHART=${metrics_dir}/{1}_gc_bias_metrics.pdf \
      S=${metrics_dir}/{1}_gc_summary_metrics.txt

#Hs metrics
ls | grep _recal.bam | sed 's/_recal.bam//' | parallel --max-args=1 \
java -Xmx10G -jar ${picard}/picard.jar CollectHsMetrics \
      I={1}_recal.bam \
      R=${db}/GRCh38_character_genome.fasta \
      O=${metrics_dir}/{1}_hs_metrics.txt \
      BAIT_INTERVALS=${db}/hg38_no_alt.interval_list \
      TARGET_INTERVALS=${db}/hg38_no_alt.interval_list \
      PER_BASE_COVERAGE=${metrics_dir}/{1}_per_base_coverage.txt
      
#VCF metrics
ls | grep _recal.bam | sed 's/_recal.bam//' | parallel --max-args=1 \
java -Xmx8G -jar ${picard}/picard.jar CollectVariantCallingMetrics \
      I={1}_final.vcf.gz  \
      O=${metrics_dir}/{1}_vcf_metrics \
      DBSNP=${db}/dbsnp_156.vcf.gz \
      TARGET_INTERVALS=${db}/abanalitica_exome.interval_list

#%PASS, %PASS (DP>=10)
#nano ${metrics_dir}/prop_var_dp10.txt
#${scripts}/calculate_proportion.sh X4493_S1_final.vcf.gz 10

# Create the CSV file and add the header
output_csv="${metrics_dir}/proportions.csv"
echo "Filename;Proportion" > "$output_csv"

# Iterate over all _final.vcf.gz files in the current directory
for vcf_file in *_final.vcf.gz; do
    if [ -f "$vcf_file" ]; then
        proportion=$(${scripts}/calculate_proportion.sh "$vcf_file" 10)
        echo "$vcf_file;$proportion" >> "$output_csv"
    fi
done

echo "Results saved in $output_csv"


