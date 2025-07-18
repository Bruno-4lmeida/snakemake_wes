rule bwa_mem_2:
    input:
        r1=config["output"] + "output_interomics/fastp/{sample}_fastp_1.fastq.gz",
        r2=config["output"] + "output_interomics/fastp/{sample}_fastp_2.fastq.gz",
    output:
        bam=config["output"] + "output_interomics/aligned/{sample}_aligned.bam",
    conda:
        "../envs/bwa.yaml"  # Certifique-se de que bwa-mem2 esteja incluído neste ambiente
    params:
        threads=13,
        ref=config["refgenome"]["ref"],
        dir=config["output"] + "output_interomics/aligned/",
        temp=config["output"] + "output_interomics/aligned/temp/",
        out=config["output"] + "output_interomics/aligned/{sample}_aligned.bam",
        sample_name="{sample}",
        read_group="'@RG\\tID:{sample}\\tSM:{sample}\\tLB:SureSelectV7\\tPL:ILLUMINA'"
    benchmark:
        config["output"] + "output_interomics/benchmarks/aligned/bwa2_{sample}_benchmark.txt",
    
    shell:
        """
        echo "##############################################" 
        echo "----------- Running BWA MEM2 ---------------" 
        echo "##############################################"

        ## Criar diretórios para saída temporária
        mkdir -p {params.dir}
        mkdir -p {params.temp}

        ## Executar bwa-mem2
        bwa-mem2 mem \
        -t {params.threads} \
        -R {params.read_group} \
        {params.ref} \
        {input.r1} {input.r2} | samtools sort --threads {params.threads} -o {params.out} -T {params.temp}
        """
