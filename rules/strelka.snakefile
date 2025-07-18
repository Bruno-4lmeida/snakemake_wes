rule strelka:
    input:
        bam=config["output"] + "output_interomics/recal_apply/{sample}_markduplicates_recal.bam",
        ref=config["refgenome"]["ref"],
        known_var=config["refgenome"]["known_var"],
        bed=config["refgenome"]["bed"],
    output:
        vcf=config["output"] + "output_interomics/vcf/call/strelka/{sample}.vcf.gz",
    conda:
        "../envs/strelka.yaml"
    params:
        out=config["output"] + "output_interomics/vcf/call/strelka/variants_{sample}",
        bed=config["refgenome"]["bed"],
        temp_dir=config["output"] + "output_interomics/temp/{sample}",
        known_var=config["refgenome"]["known_var"],
    threads: 56
    shell:
         """
        # Configurar o workflow do Strelka
        configureStrelkaGermlineWorkflow.py \
            --bam {input.bam} \
            --referenceFasta {input.ref} \
            --runDir {params.temp_dir} \
            --callRegions {input.bed}
        
        {params.temp_dir}/runWorkflow.py -m local -j {threads}

        mv {params.temp_dir}/results/variants/genome.S1.vcf.gz {params.out}
        mv {params.temp_dir}/results/variants/genome.S1.vcf.gz.tbi {params.out}.tbi

        mv {params.temp_dir}/results/variants/variants.vcf.gz {output.vcf}
        mv {params.temp_dir}/results/variants/variants.vcf.gz.tbi {output.vcf}.tbi
        
        """
