rule mutect_somatic:
    input:
        tumor_bam = lambda wildcards: f"{config['output']}output_interomics/recal_apply/{wildcards.tumor}_markduplicates_recal.bam",
        normal_bam = lambda wildcards: f"{config['output']}output_interomics/recal_apply/{wildcards.normal}_markduplicates_recal.bam"
    output:
        vcf = config["output"] + "output_interomics/mutect/{tumor}_vs_{normal}_filtered.vcf.gz"
    params:
        ref = config["refgenome"]["ref"],
        resources = config["refgenome"]["known_var"],
        intervals = config["refgenome"]["interval"],
    threads: 4
    shell:
        """
        {GATK}/gatk Mutect2 \
            -R {params.ref} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -tumor {wildcards.tumor} \
            -normal {wildcards.normal} \
            --germline-resource {params.resources} \
            --intervals {params.intervals} \
            -O {output.vcf}


        ls | grep _recal.bam | sed 's/_.*//' | parallel --max-args=1 \
        ${GATK}/gatk GetPileupSummaries \
        -I {1}_recal.bam  \
        -V ${GATK_resources}/af-only-gnomad.hg38.vcf \
        -L ${GATK_resources}/hg38.bed \
        -O ../variant_calling/{1}_pileups.table

        cat ${tmpgenome}/somatic_wgs.txt | awk -F" " 'NR>=1 {print $1}' | parallel  -n 2 -j 4 \
        ${GATK}/gatk CalculateContamination \
        -I ../variant_calling/{1}_pileups.table \
        -matched ../variant_calling/{2}_pileups.table \
        -O ../variant_calling/{2}_contamination.table \
        --tumor-segmentation ../variant_calling/{2}_segments.table

        for i in $(ls | grep _normal_tumor.vcf.gz | sed 's/_.*//' | uniq); do

        #replacing NaN for 0.0 to the contamination table file when needed
        cat $i"_contamination.table" | sed 's/NaN/0.0/' > $i"_edited.contamination.table"

        ${GATK}/gatk FilterMutectCalls \
        -R ${db}/GRCh38_genome.fa \
        -V $i"_normal_tumor.vcf.gz" \
        --contamination-table $i"_edited.contamination.table" \
        --tumor-segmentation $i"_sorted.segments.table" \
        -O $i"_filtered.vcf.gz"

        #### CALCULATE VAF####
        bcftools +fill-tags $i"_filtered.vcf.gz" \
        -O z \
        -o $i"_filtered_VAF.vcf.gz" \
        -- -t VAF

        mv $i"_filtered_VAF.vcf.gz" $i"_filtered.vcf.gz" 


        """