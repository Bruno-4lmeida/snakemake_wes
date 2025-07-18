rule process_hs_metrics:
    input:
        file=expand(config["output"] + "output_interomics/metrics/{sample}_hs_metrics.txt", sample=SAMPLES),
    output:
        config["output"] + "output_interomics/metrics/summary/consolidated_hs_metrics.csv"
    params:
        summ_dir=config["output"] + "output_interomics/metrics/summary/",
        met_dir=config["output"] + "output_interomics/metrics/"

    conda:
        "../envs/r_metrics.yaml"
    shell:
        """
        mkdir -p {params.summ_dir}
        mkdir -p {params.met_dir}

        Rscript --vanilla /home/operator/InterOmics/workflows/interomics_snakemake/scripts/process_hs_metrics.R "{params.met_dir}" "{output}"
        """
