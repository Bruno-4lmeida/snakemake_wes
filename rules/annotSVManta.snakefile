rule annotSVManta:
    input:
        vcf_sv=config["output"] + "output_interomics/manta/{sample}/{sample}_fixed.vcf",
    output:
        vcf_sv=config["output"] + "output_interomics/annotsv/manta/{sample}_annotsv_manta.tsv",
    params:
        output_dir=config["output"] + "output_interomics/annotsv/manta/",
        db=config["refgenome"]["annotsvdb"],  # Look this path and change when it needs
    threads: 30
    conda:
        "../envs/annotsv.yaml"
    shell:
        """
        echo "##############################################" 
        echo "-----------    Running AnnotSV    ------------" 
        echo "##############################################"           

        ## SV annotation
        AnnotSV \
            -SVinputFile {input.vcf_sv} \
            -outputFile {output.vcf_sv} \
            -annotationsDir {params.db}/AnnotSV_annotations  \
            -genomeBuild GRCh38  

        """
