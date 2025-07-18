rule strelka:
    input:
        normal_bam=config["output"] + "output_interomics/recal_apply/{sample}_N_recal.bam",
        tumor_bam=onfig["output"] + "output_interomics/recal_apply/{sample}_T_recal.bam",
        ref=config["refgenome"]["ref"]
    output:
        vcf=config["output"] + "output_interomics/vcf/call/strelka/{sample}.vcf.gz",
    conda:
        "../envs/strelka.yaml"
    params:
        out=config["output"] + "output_interomics/vcf/call/strelka/variants_{sample}",
        bed=config["refgenome"]["bed"],
        temp_dir_manta=config["output"] + "output_interomics/temp_manta/{sample}_possibleindels",
        temp_dir=config["output"] + "output_interomics/temp/{sample}",
    threads: 56 
    shell:
         """
        #manta
        configManta.py \
        --normalBam {input.normal_bam} \
        --tumorBam {input.tumor_bam} \
        --referenceFasta {input.ref} \
        --runDir {params.temp_dir_manta} \
        --exome

        {params.temp_dir_manta}/runWorkflow.py -m local -j {threads}

        # Configurar o workflow do Strelka
        configureStrelkaSomaticWorkflow.py \
        --normalBam {input.normal_bam} \
        --tumorBam  {input.tumor_bam} \
        --referenceFasta {input.ref}  \
        --indelCandidates {params.temp_dir_manta}/results/variants/candidateSmallIndels.vcf.gz \
        --exome \
        --runDir {params.temp_dir}
        
        {params.temp_dir}/runWorkflow.py -m local -j {threads}

        mv {params.temp_dir}/results/variants/genome.S1.vcf.gz {params.out}
        mv {params.temp_dir}/results/variants/genome.S1.vcf.gz.tbi {params.out}.tbi

        mv {params.temp_dir}/results/variants/variants.vcf.gz {output.vcf}
        mv {params.temp_dir}/results/variants/variants.vcf.gz.tbi {output.vcf}.tbi
        
        """
