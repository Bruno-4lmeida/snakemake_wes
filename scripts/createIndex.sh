conda activate picard

# CREATE REF.dict
PicardCommandLine CreateSequenceDictionary R=GRCh38_numeric_genome.fasta O=GRCh38_numeric_genome.dict
PicardCommandLine CreateSequenceDictionary R=GRCh38_character_genome.fasta O=GRCh38_character_genome.dict

# create file.fai
samtools faidx GRCh38_numeric_genome.fasta
samtools faidx GRCh38_character_genome.fasta

# BWA INDEX
bwa index GRCh38_numeric_genome.fasta
bwa index GRCh38_character_genome.fasta

## CREATE SAMPLES TEST
zcat COL53_1.fastq.gz | head -n 4000 > /home/operator/InterOmics/workflows/interomics_snakemake/samples/COL53_1_subset_reads_1.fastq
pigz /home/operator/InterOmics/workflows/interomics_snakemake/samples/COL53_1_subset_reads_1.fastq

zcat COL53_2.fastq.gz | head -n 4000 > /home/operator/InterOmics/workflows/interomics_snakemake/samples/COL53_1_subset_reads_2.fastq
pigz /home/operator/InterOmics/workflows/interomics_snakemake/samples/COL53_1_subset_reads_2.fastq
