rule sansa:
    input:
        vcf_cnv=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.vcf.gz",
        vcf_sv=config["output"] + "output_interomics/delly/svs/{sample}_sv.vcf.gz",
        gnomad=config["refgenome"]["gnomad_sv"],
        gtf=config["refgenome"]["gtf"],
    output:
        cnv_anno=config["output"] + "output_interomics/sansa/cnvs/{sample}_cnv_anno.bcf",
        cnv_query=config["output"] + "output_interomics/sansa/cnvs/{sample}_cnv_query.tsv.gz",
        sv_anno=config["output"] + "output_interomics/sansa/svs/{sample}_SV_anno.bcf",
        sv_query=config["output"] + "output_interomics/sansa/svs/{sample}_SV_query.tsv.gz",

    conda:
        "../envs/sansa.yaml"
    threads: 3
    shell:
        """
        echo "##############################################"
        echo "------  Running SANSA  SV annotation   -------"
        echo "##############################################"  

        ## -S report all matches and -m keep unmatcheds
        sansa annotate -m -c -s all -g {input.gtf} -d {input.gnomad} {input.vcf_cnv} -a {output.cnv_anno} -o {output.cnv_query}
        sansa annotate -m -c -s all -g {input.gtf} -d {input.gnomad} {input.vcf_sv} -a {output.sv_anno} -o {output.sv_query}

        """
