rule recal_table:
    input:
        bam=config["output"] + "output_interomics/marked/{sample}_markduplicates.bam",
        ref=config["refgenome"]["ref"],
        dict=config["refgenome"]["dict"],
        known_var=config["refgenome"]["known_var"],  # optional known sites - single or a list
    output:
        recal_table=config["output"] + "output_interomics/recal_table/{sample}_recal.grp",
    log:
        config["output"] + "output_interomics/recal_table/log/{sample}.log",
    conda:
        "../envs/gatk.yaml"
    threads: config["threads"]["recal_table"],
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        #mem_mb=config["ram_memory"]["recal_table"],
        mem_mb=25000,
    benchmark: config["output"] + "output_interomics/benchmarks/recal_table/recal_table_{sample}_benchmark.txt",

    shell:
        """
        echo "##############################################" 
        echo "-----   Running GATK BaseRecalibrator  -------" 
        echo "##############################################"

        gatk BaseRecalibrator \
        -I {input.bam} \
        -R {input.ref} \
        --known-sites {input.known_var} \
        -O {output.recal_table}
        """

