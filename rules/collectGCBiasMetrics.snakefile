rule CollectGcBiasMetrics:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
    output:
        gc_bias=config["output"] + "output_interomics/metrics/{sample}_gc_bias_metrics.txt",
        gc_pdf=config["output"] + "output_interomics/metrics/{sample}_gc_bias_metrics.pdf",
        gc_summ=config["output"] + "output_interomics/metrics/{sample}_gc_summary_metrics.txt",
    conda:
        "../envs/picard.yaml"
    threads: 2
    params:
        ref=config["refgenome"]["ref"],
        

    shell:
        """
        echo "##############################################" 
        echo "-----  Picard CollectGCBiasMetrics     -------" 
        echo "##############################################"
        picard CollectGcBiasMetrics \
            I={input.bam} \
            R={params.ref} \
            O={output.gc_bias} \
            CHART={output.gc_pdf} \
            S={output.gc_summ}
        sync # forcar o snakemake a sincronizar os arquivos.
        """
