rule trim_galore_single:
    input:
        config["fastp"]["dir"] + "{sample}_fastp.fastq"
    output:
        config["trim-galore"]["dir"] + "{sample}_trimmed.fastq"
    params:
        threads=config["threads"],
        dir=directory(config["trim-galore"]["dir"]),
    log:
        config["trim-galore"]["dir"] + "log/{sample}.log",
    conda:
        "../envs/trim_galore.yaml"
    shell:
        """
        echo "##############################################"
        echo "------------    Running trim-galore  ---------"
        echo "##############################################"

        trim-galore {params.dir} {input} 2> {log}
        """




