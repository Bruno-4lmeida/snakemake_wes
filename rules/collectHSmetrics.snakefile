rule HSMetrics:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
    output:
        out_hs=config["output"] + "output_interomics/metrics/{sample}_hs_metrics.txt",
    conda:
        "../envs/picard.yaml"
    threads: 2
    params:
        dir=config["output"] + "output_interomics/metrics/{sample}_per_base_coverage.txt.gz",
        ref=config["refgenome"]["ref"],
        intervals=config["refgenome"]["interval"],
    shell:
        """
        echo "##############################################" 
        echo "-----  PicardCommandLine CollectHsMetrics ----" 
        echo "##############################################"

        picard CollectHsMetrics \
            I={input.bam} \
            R={params.ref} \
            O={output.out_hs} \
            BAIT_INTERVALS={params.intervals} \
            TARGET_INTERVALS={params.intervals} \
            PER_BASE_COVERAGE=/dev/stdout | gzip > {params.dir}
        """
