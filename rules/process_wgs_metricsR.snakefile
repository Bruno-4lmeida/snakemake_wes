rule process_metrics:
    input:
        files=expand(config["output"] + "output_interomics/metrics/{sample}_collect_wgs_metrics.txt", sample=SAMPLES),
    output:
       config["output"] + "output_interomics/metrics/summary/consolidated_wgs_metrics.csv"
    params:
        sum_dir=config["output"] + "output_interomics/metrics/summary/",
        met_dir=config["output"] + "output_interomics/metrics/"


    conda:
        "../envs/r_metrics.yaml"
    shell:
        """
        mkdir -p {params.sum_dir}
        mkdir -p {params.met_dir}


        Rscript --vanilla /home/operator/InterOmics/workflows/interomics_snakemake/scripts/process_wgs_metrics.R "{params.met_dir}" "{output}"
        """