rule delly_genotyping:
    input:
        cnv_bcf=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.bcf",
        sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.bcf",
    output:
        geno_cnv_bcf=config["output"] + "output_interomics/delly/cnvs/{sample}_geno_cnv.bcf",
        geno_sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_geno_sv.bcf",

    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    conda:
        "../envs/delly.yaml"
    threads: 5
    shell:
        """
        echo "##############################################"
        echo "--------  Running genotyping Delly -----------"
        echo "##############################################"  

        # Genotipar CNVs
        delly cnv -u -v {input.sv_bcf} -g {params.reference} -m {params.mappability} -o {output.geno_cnv_bcf} {input.cnv_bcf}

        # Genotipar SVs
        delly call -u -v {input.sv_bcf} -g {params.reference} -o {output.geno_sv_bcf} {input.sv_bcf}

        """
