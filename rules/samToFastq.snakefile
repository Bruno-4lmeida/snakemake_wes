rule sam_to_bam:
    input:
        config["picard"]["dir"] + "{sample}_picard.bam"
    output:
        config["samToFastq"]["dir"] + "{sample}_clipped.fastq"
    
#For single end
    shell:
        """
        echo "##############################################" 
        echo "---------- PICARD SamToFastq    --------------" 
        echo "##############################################"

        PicardCommandLine SamToFastq  \
        --INPUT {input} \
        --FASTQ {output} \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2 \
        --INTERLEAVE false \
        --INCLUDE_NON_PF_READS true
        """