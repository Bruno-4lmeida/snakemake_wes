rule cnvkit_call:
    input:
        bam="config["output"] + "output_vkaryo/marked/{sample}_markduplicates.bam",
        reference=config["refgenome"]["ref"],,
        targets="reference/targets.bed",
        antitargets="reference/antitargets.bed"
    output:
        cnr=config["output"] + "output_vkaryo/cnvkit/{sample}.cnr",
        cns=config["output"] + "output_vkaryo/cnvkit/{sample}.cns"
    params:
        access="reference/access.bed",
        output_dir=config["output"] + "output_vkaryo/cnvkit/"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        """
        # Criação de uma referência a partir dos arquivos BED de alvos e antialvos
        cnvkit.py reference {input.bam} -t {input.targets} -a {input.antitargets} -f {input.reference} -o {params.output_dir}/reference.cnn

        # Correção de cobertura
        cnvkit.py fix {params.output_dir}/{wildcards.sample}.targetcoverage.cnn {params.output_dir}/{wildcards.sample}.antitargetcoverage.cnn {params.output_dir}/reference.cnn -o {output.cnr}

        # Segmentação
        cnvkit.py segment {output.cnr} -o {output.cns}

        # Chamada de CNVs
        cnvkit.py call {output.cns} -o {params.output_dir}/{wildcards.sample}_call.cns
        """
