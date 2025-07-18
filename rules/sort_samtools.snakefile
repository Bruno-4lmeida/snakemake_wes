rule sort_samtools:
    input:
        config["bwa"]["dir"] + "{sample}_aligned.bam"
    output:
        config["samtools"]["dir"] + "{sample}_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo "##############################################" 
        echo "--------    Running SAMTOOLS Sort   ----------" 
        echo "##############################################"
        samtools sort {input} -T {wildcards.sample} -o {output}
        """
