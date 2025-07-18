rule delly_cnv:
    input:
        bam=config["input"] + "{sample}_recal.bam",
        sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.bcf",
    output:
        cnv_bcf=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.bcf",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    threads: 5
    conda:
        "../envs/delly.yaml"
    shell:
        """
        echo "##############################################"
        echo "--------    Running DELLY CNV    ------------"
        echo "##############################################"  

        delly cnv -g {params.reference} -m {params.mappability} -o {output.cnv_bcf} -l {input.sv_bcf} {input.bam}
        """
