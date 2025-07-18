rule fastq_to_bam:
    input:
        config["fastp"]["dir"] + "{sample}_fastp.fastq"
    output:
        config["picard"]["dir"] + "{sample}_picard.bam"

    shell:
        """
        PicardCommandLine FastqToSam \
            FASTQ={input} \
            OUTPUT={output} \
            SAMPLE_NAME={wildcards.sample} \
            TMP_DIR=tmp
        """
