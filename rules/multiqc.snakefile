# multiqc - merge fastqc reports
# *********************************************************************
rule multiqc:
    input:
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp.html", sample=SAMPLES),
        expand(config["output"] + "output_interomics/fastp/{sample}_fastp.json", sample=SAMPLES),
    output:
        dir=directory(config["output"] + "output_interomics/multiqc/"),
    params:
        config["output"] + "output_interomics/fastp/"
    log:
        config["output"] + "output_interomics/multiqc/log/multiqc.log",
    conda:
        "../envs/multiqc.yaml",
    threads: config["threads"]["others"],

    shell:
        """
        echo "##############################################"
        echo "-----------    Running MultiQC    ------------"
        echo "##############################################"
        
        multiqc --force --outdir {output.dir} {input} 2> {log}
        """