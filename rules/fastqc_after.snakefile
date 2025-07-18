rule fastqc_after_single:
    input:
        "data/fastp/{sample}_fastp.fastq"
    output:
        "data/fastqc/{sample}_afterfastqc.html"
    shell:
        "fastqc -o data/fastqc {input}"