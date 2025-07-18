import math
from snakemake.io import glob_wildcards

(SAMPLES,) = glob_wildcards(config["input"] + "{sample}_1.fastq.gz")
num_amostras = len(SAMPLES)
num_cores = workflow.cores
rule fastp:
    input:
        r1=config["input"] + "{sample}_1.fastq.gz",
        r2=config["input"] + "{sample}_2.fastq.gz"
    output:
        r1=config["output"] + "output_interomics/fastp/{sample}_fastp_1.fastq.gz",
        r2=config["output"] + "output_interomics/fastp/{sample}_fastp_2.fastq.gz",
        html=config["output"] + "output_interomics/fastp/{sample}_fastp.html",
        json=config["output"] + "output_interomics/fastp/{sample}_fastp.json",
    log:
        config["output"] + "output_interomics/fastp/log/{sample}.log"
    threads: math.ceil(num_cores/ num_amostras),
    params:
        threads=config["threads"],
        dir=directory(config["output"] + "output_interomics/fastp/")
    conda:
        "../envs/fastp.yaml"
    benchmark: config["output"] + "output_interomics/benchmarks/fastp/fastp_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################"
        echo "------------    Running Fastp    ------------"
        echo "##############################################"

        fastp --thread {threads} -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --html {output.html} --json {output.json}
        """