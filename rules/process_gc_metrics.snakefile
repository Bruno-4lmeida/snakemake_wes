rule process_gc_metrics:
    input:
        file=expand(config["output"] + "output_interomics/metrics/{sample}_gc_summary_metrics.txt", sample=SAMPLES),
    output:
        config["output"] + "output_interomics/metrics/summary/consolidated_gc_metrics.csv"
    params:
        config["output"] + "output_interomics/metrics/summary/",
    conda:
        "../envs/r_metrics.yaml"
    shell:
        """
        mkdir -p {params}
        Rscript --vanilla /home/operator/InterOmics/workflows/interomics_snakemake/scripts/process_gc_metrics.R "{params}" "{output}"
        """
