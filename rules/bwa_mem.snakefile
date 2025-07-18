rule bwa_mem:
    input:
       r1=config["output"] + "output_interomics/fastp/{sample}_fastp_1.fastq.gz",
       r2=config["output"] + "output_interomics/fastp/{sample}_fastp_2.fastq.gz",
    output:
        bam=config["output"] + "output_interomics/aligned/{sample}_aligned.bam",
    conda:
        "../envs/bwa.yaml"
    threads: config["threads"]["bwa"]
    params:
        threads= config["threads"]["bwa"],
        ref=config["refgenome"]["ref"],
        dir=config["output"] + "output_interomics/aligned/",
        temp=config["output"] + "output_interomics/aligned/temp/",
        out=config["output"] + "output_interomics/aligned/{sample}_aligned.bam",
        sample_name="{sample}",
        read_group="'@RG\\tID:{sample}\\tSM:{sample}\\tLB:SureSelectV7\\tPL:ILLUMINA'"

    benchmark: config["output"] + "output_interomics/benchmarks/aligned/bwa_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################" 
        echo "-----------    Running BWA MEM    ------------" 
        echo "##############################################"
        ## Create read_group
        #fastq="{input.r1}"
        ##header=$(zcat -f $fastq | head -n 1)
        #header=$(head -n 1 $fastq)
        #read_id=${{header#@}}
        #sample_name=${{read_id}}  # Usando o ID da leitura como nome da amostra
        #library_name=${{read_id}}  # Usando o ID da leitura como nome da biblioteca
        #platform="ILLUMINA"
        ##platform_unit1="${{read_id}}.1"
        ##read_group="@RG\\tID:${{read_id}}\\tSM:${{sample_name}}\\tPL:${{platform}}\\tLB:${{library_name}}\\tPU:${{platform_unit}}"
        #read_group="@RG\\tID:${{read_id}}\\tSM:${{sample_name}}\\tPL:${{platform}}"

        ## running bwa
        mkdir -p {params.dir}
        mkdir -p {params.temp}

        #cd {params.dir}

        #bwa mem -t {threads} {params.ref} -R "$read_group" {input.r1} {input.r2} -o {params.out}

        bwa mem \
        -M \
        -t {params.threads} \
        -v 3 \
        -R {params.read_group} \
        {params.ref} \
        {input.r1} {input.r2} | samtools sort --threads {params.threads} -o {params.out} -T {params.temp}

        """