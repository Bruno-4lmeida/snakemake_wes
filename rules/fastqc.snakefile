rule fastqc_before:
    input:
        r1=config["input"]["samples"] + "{sample}_1.fastq.gz",
        r2=config["input"]["samples"] + "{sample}_2.fastq.gz",
    params:
        out_dir=config["output"] + config["fastqc"]
    output:
        dir=directory(config["fastqc"]["dir"]),
        fastqc_html_1 = expand(config["output"] + config["fastqc"] + "{sample}_1_fastqc.html", sample=SAMPLES),
        fastqc_zip_1 = expand(config["fastqc"]["dir"] + "{sample}_1_fastqc.zip", sample=SAMPLES),
        fastqc_html_2 = expand(config["fastqc"]["dir"] + "{sample}_2_fastqc.html", sample=SAMPLES),
        fastqc_zip_2 = expand(config["fastqc"]["dir"] + "{sample}_2_fastqc.zip", sample=SAMPLES),  
    conda:
        "../envs/fastqc.yaml"
    threads: config["threads"]["others"],
    shell:
        """
        echo "##############################################" 
        echo "------------    Running FastQC    ------------" 
        echo "##############################################"
        mkdir -p {params.out_dir}
        fastqc -o {params.out_dir} {input.r1} {input.r2} 
        """




    