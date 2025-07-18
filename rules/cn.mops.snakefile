rule cnmops:
    input:
        expand(config["input"] + "{sample}_markduplicates_recal.bam", sample=SAMPLES) # tem que trocar aqui input p output
    output:
        cnv_output=config["output"] + "output_interomics/cnvs/cnmops/resultscnmops_cnvs.csv"
    params:
        dir=config["input"] + "output_interomics/recal_apply/", # trocar aqui tbm input p output
        reference=config["refgenome"]["ref"],
        bed=config["refgenome"]["bed"],
    threads: 50
    conda:
        "../envs/r_metrics.yaml"  
    shell:
        """
        Rscript  --vanilla /home/operator/InterOmics/workflows/interomics_snakemake/scripts/cn.mops.R "{params.dir}" "{output.cnv_output}"  "{threads}" "{params.bed}"

        """
