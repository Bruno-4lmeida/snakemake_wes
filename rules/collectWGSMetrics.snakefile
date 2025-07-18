rule wgsMetrics:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
    output:
        wgs_metrics=config["output"] + "output_interomics/metrics/{sample}_collect_wgs_metrics.txt",
    conda:
        "../envs/picard.yaml"
    threads: 2
    params:
        ref=config["refgenome"]["ref"],
        intervals=config["refgenome"]["interval"],
    shell:
        """
        echo "##############################################" 
        echo "-----  METRICS  -------" 
        echo "##############################################"

        picard CollectWgsMetrics \
            I={input.bam} \
            O={output.wgs_metrics} \
            R={params.ref} \
            INTERVALS={params.intervals} \
            MINIMUM_MAPPING_QUALITY=20 \
            MINIMUM_BASE_QUALITY=20 \
            COVERAGE_CAP=250 

        """
