rule markDuplicates:
    input:
        config["output"] + "output_interomics/aligned/{sample}_aligned.bam",
    output:
        bam=config["output"] + "output_interomics/marked/{sample}_markduplicates.bam",
        txt=config["output"] + "output_interomics/marked/{sample}_markduplicates_metrics.txt",
    conda:
        "../envs/picard.yaml"
    threads: config["threads"]["others"],
    resources:
        #mem_mb=config["ram_memory"]["markduplicates"],
        mem_mb=25000,
    benchmark: config["output"] + "output_interomics/benchmarks/mark_duplicates/marked_{sample}_benchmark.txt",


    shell:
        """
        echo "##############################################" 
        echo "-----   Running Picard MarkDuplicates   -----" 
        echo "##############################################"

        picard MarkDuplicates \
        I={input} \
        O={output.bam}\
        M={output.txt} \
        REMOVE_DUPLICATES=true \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ADD_PG_TAG_TO_READS=false \
        PROGRAM_RECORD_ID=null \
        ASSUME_SORT_ORDER=coordinate
        """