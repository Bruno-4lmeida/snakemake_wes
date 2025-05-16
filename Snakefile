# MAIN SNAKEFILE

configfile: "config.yaml"
# To get the samples names form sample dir
(SAMPLES,) = glob_wildcards(config["input"] + "{sample}_1.fastq.gz")
#(SAMPLES,) = glob_wildcards("/home/operator/Downloads/ab_analitica/output_interomics/recal_apply/{sample}_markduplicates_recal.bam")
#SAMPLES=['gDNA205', 'gDNA254', 'COL53M', 'DNA475_S10', 'DNA405_S4', 'MT23_S1', 'COL40M', 'X22_2310_S9', 'DNA491_S8', 'COL40', 'MT23M_S2',  'DNA485_S7', 'COL40P', 'COL53P', 'gDNA93', 'DNA494_S5', 'X23_3080_S7', 'DNA490_S4', 'COL56_S10', 'DNA406_S5', 'DNA412_S6', 'MT23P_S3', 'DNA500_S6', 'COL53', 'X23_3083_S8', 'DNA492_S9']
print(SAMPLES)

print(f"Sample names: {', '.join(SAMPLES)}")

if config["use_bwa"]:
    include: "rules/bwa_mem.snakefile"
else: 
    include: "rules/dragmap.snakefile"

if config["germline"]:
    if config["use_gatk"]:
        include: "rules/haplotypeCaller.snakefile"
        include: "rules/genotyping_gatk.snakefile"

    if config["use_strelka"]:  # Se for para usar o Strelka
        include: "rules/strelka.snakefile"

    if config["use_deepvariant"]: 
        import subprocess
        try:
            # Verifica se o Docker está instalado
            subprocess.run(['singularity', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            singularity_installed = True
        except subprocess.CalledProcessError:
            singularity_installed = False

        if singularity_installed:
            print("Singularity is installed.")
            include: "rules/deepvariant_singularity.snakefile"
        else:
            print(
                """
                    ERROR: Singularity is not installed or is not accessible!
                    Install singularity and try again or use another
                    caller!
                """)
else:
    include: "rules/strelka_somatic.snakefile"

# To import the rules for each tool
include: "rules/multiqc.snakefile"
include: "rules/fastp.snakefile"
include: "rules/markduplicates.snakefile"
include: "rules/recal_table.snakefile"
include: "rules/recal_apply.snakefile"
include: "rules/selectVar.snakefile"
include: "rules/filterVar.snakefile"
include: "rules/concat_bcftools.snakefile"
include: "rules/annovar.snakefile"
include: "rules/collectGCBiasMetrics.snakefile"
include: "rules/collectHSmetrics.snakefile"
include: "rules/collectWGSMetrics.snakefile"
include: "rules/collectVariantCallingMetrics.snakefile"
include: "rules/collect_proportions.snakefile"
include: "rules/process_wgs_metricsR.snakefile"
include: "rules/process_gc_metrics.snakefile"
include: "rules/process_hs_metrics.snakefile"


# rule all for all output files
rule all:
    input:
        config["output"] + "output_interomics/multiqc/", # output from multiqc
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp_1.fastq.gz", sample=SAMPLES), # output from fastp
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp_2.fastq.gz", sample=SAMPLES), # output from fastp
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp.html", sample=SAMPLES), # output from fastp
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp.json", sample=SAMPLES),   # output from fastp
        ### ALIGNED; MARKED; RECAL_APPLY
        expand(config["output"] + "output_interomics/aligned/{sample}_aligned.bam", sample=SAMPLES), # output from mappers
        expand(config["output"] + "output_interomics/marked/{sample}_markduplicates.bam", sample=SAMPLES), # output from markduplicates        
        expand(config["output"] + "output_interomics/marked/{sample}_markduplicates_metrics.txt", sample=SAMPLES), # output from markduplicates
        expand(config["output"] + "output_interomics/recal_table/{sample}_recal.grp", sample=SAMPLES), # output from baserecalibration
        expand(config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam", sample=SAMPLES), #output base recalibration (rule recal_apply)
#        ### VARIANT CALLER ###
        (expand(config["output"] + "output_interomics/vcf/call/haplotype_caller/{sample}.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]),       
        (expand(config["output"] + "output_interomics/vcf/call/strelka/{sample}.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]),
        (expand(config["output"] + "output_interomics/vcf/call/deep_variant/{sample}.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]),
        ### SELECT ###
        (expand(config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_snps.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]),#output selectVar
        (expand(config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_snps.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]),#output selectVar
        (expand(config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_snps.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]),#output selectVar
        (expand(config["output"] + "output_interomics/vcf/selected/haplotype_caller/{sample}_raw_indels.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]),#output selectVar
        (expand(config["output"] + "output_interomics/vcf/selected/strelka/{sample}_raw_indels.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]),#output selectVar
        (expand(config["output"] + "output_interomics/vcf/selected/deep_variant/{sample}_raw_indels.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]),#output selectVar
        ## FILTERING ###
        (expand(config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_snps.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]),#output filterVar
        (expand(config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_snps.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]),#output filterVar
        (expand(config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_snps.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]),#output filterVar
        (expand(config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered_indels.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]),#output filterVar
        (expand(config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered_indels.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]),#output filterVar
        (expand(config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered_indels.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]),#output filterVar
        ### CONCAT ###
        (expand(config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_filtered.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]), #output concat
        (expand(config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_filtered.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]), #output concat
        (expand(config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_filtered.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]), #output concat
        (expand(config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz", sample=SAMPLES) if config["use_gatk"] == True else[]), #output concat
        (expand(config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz", sample=SAMPLES) if config["use_strelka"] == True else[]), #output concat
        (expand(config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz", sample=SAMPLES) if config["use_deepvariant"] == True else[]), #output concat
        ### ANNOTATION ###
        (expand(config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annotation_updated.hg38_multianno.csv", sample=SAMPLES) if config["use_gatk"] == True else[]), #output annovar
        (expand(config["output"] + "output_interomics/annovar/strelka/{sample}_annotation_updated.hg38_multianno.csv", sample=SAMPLES) if config["use_strelka"] == True else[]), #output annovar
        (expand(config["output"] + "output_interomics/annovar/deep_variant/{sample}_annotation_updated.hg38_multianno.csv", sample=SAMPLES) if config["use_deepvariant"] == True else[]), #output annovar
        ### METRICS ###
        expand(config["output"] + "output_interomics/metrics/{sample}_collect_wgs_metrics.txt", sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/{sample}_gc_summary_metrics.txt",  sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/{sample}_gc_bias_metrics.pdf", sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/{sample}_hs_metrics.txt", sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/summary/proportions.csv", sample=SAMPLES),
        ### VCF METRICS
        (expand(config["output"] + "output_interomics/metrics/haplotype_caller/{sample}_vcf_metrics.variant_calling_detail_metrics", sample=SAMPLES) if config["use_gatk"] == True else[]),
        (expand(config["output"] + "output_interomics/metrics/strelka/{sample}_vcf_metrics.variant_calling_detail_metrics", sample=SAMPLES) if config["use_strelka"] == True else[]),
        (expand(config["output"] + "output_interomics/metrics/deep_variant/{sample}_vcf_metrics.variant_calling_detail_metrics", sample=SAMPLES) if config["use_deepvariant"] == True else[]),

        ### PROCESS METRICS ###
        expand(config["output"] + "output_interomics/metrics/summary/consolidated_wgs_metrics.csv", sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/summary/consolidated_gc_metrics.csv", sample=SAMPLES),
        expand(config["output"] + "output_interomics/metrics/summary/consolidated_hs_metrics.csv", sample=SAMPLES),
