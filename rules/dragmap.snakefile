rule dragmap:
    input:
       r1=config["output"] + "output_interomics/fastp/{sample}_fastp_1.fastq.gz",
       r2=config["output"] + "output_interomics/fastp/{sample}_fastp_2.fastq.gz",
    output:
        bam=config["output"] + "output_interomics/aligned/{sample}_aligned.bam"
    conda:
        "../envs/dragmap.yaml"
    params:
        ref=config["refgenome"]["ref_dir"],
        temp=config["output"] + "output_interomics/aligned/temp/",
    threads: 13

    benchmark: config["output"] + "output_interomics/benchmarks/aligned/benchmarks/dragmap_{sample}_benchmark.txt"
    conda:
        "../envs/dragmap.yaml"
    shell:
        """
        echo "##############################################" 
        echo "-----------    Running dragmap    ------------" 
        echo "##############################################"

       dragen-os -r {params.ref} -1 {input.r1} -2 {input.r2} | samtools sort --threads {threads} -o {output.bam} -T {params.temp}

        

        """
