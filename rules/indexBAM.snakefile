rule indexbam:
    input:
        config["bwa"]["dir"] + "{sample}_sorted.bam" # Modificar para o ultimo bam... Sao muitos bam files
    output:
        config["bwa"]["dir"] + "{sample}_sorted.bai" # Modificar para o ultimo bam...
    
#For single end
    shell:
        """
        echo "##############################################" 
        echo "-----    Running PICARD BAM Index    ---------" 
        echo "##############################################"

        PicardCommandLine BuildBamIndex  --INPUT {input} --OUTPUT {output}
        """