rule delly_sv:
    input:
        bam=config["input"] + "{sample}_recal.bam",
    output:
        sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.bcf",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    threads: 5
    conda:
        "../envs/delly.yaml"
    shell:
        """   
        echo "##############################################"
        echo "-------    Running   DELLY SV      ----------"
        echo "##############################################"  

        delly call -g {params.reference} -o {output.sv_bcf} {input.bam}
        
        """
