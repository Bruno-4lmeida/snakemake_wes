rule delly_classify:
    input:
        cnv_bcf=config["output"] + "output_interomics/delly/cnvs/{sample}cnv.bcf",
        sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.bcf",
    output:
        cnv_bcf=config["output"] + "output_vkaryo/delly/filtered/{sample}_cnv.bcf",
        sv_bcf=config["output"] + "output_vkaryo/delly/filtered/{sample}_sv.bcf",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    conda:
        "../envs/delly.yaml"
    threads: 5
    shell:
        """
        echo "##############################################"
        echo "--------   Running Classify Delly  ------------"
        echo "##############################################"  

        # Filtrar CNVs
        delly classify -f germline -o {output.cnv_bcf} {output.cnv_bcf}

        # Filtrar SVs
        delly classify -f germline -o {output.sv_bcf} {output.sv_bcf}
        """
