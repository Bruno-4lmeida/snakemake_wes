rule annotSV:
    input:
        vcf_cnv=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.vcf.gz",
        vcf_sv=config["output"] + "output_interomics/delly/svs/{sample}_sv.vcf.gz",
    output:
        vcf_cnv=config["output"] + "output_interomics/annotsv/delly/{sample}_annotsv_cnv.tsv",
        vcf_sv=config["output"] + "output_interomics/annotsv/delly/{sample}_annotsv_SVs.tsv",
    params:
        output_dir=config["output"] + "output_interomics/annotsv/",
        db=config["refgenome"]["annotsvdb"],  # Look this path and change when it needs
    threads: 30
    conda:
        "../envs/annotsv.yaml"
    shell:
        """
        echo "##############################################" 
        echo "-----------    Running AnnotSV    ------------" 
        echo "##############################################"           

        ## CNV annotation
        AnnotSV \
            -SVinputFile {input.vcf_cnv} \
            -outputFile {output.vcf_cnv} \
            -annotationsDir {params.db}/AnnotSV_annotations  \
            -genomeBuild GRCh38

        ## SV annotation
        AnnotSV \
            -SVinputFile {input.vcf_sv} \
            -outputFile {output.vcf_sv} \
            -annotationsDir {params.db}/AnnotSV_annotations  \
            -genomeBuild GRCh38  

        """
