rule annotate_with_vep:
    input:
        cnv=config["output"] + "output_interomics/sansa/cnvs/{sample}_cnv_anno.vcf.gz",
        svs=config["output"] + "output_interomics/sansa/svs/{sample}_SV_anno.vcf.gz",
    output:
        cnv_annotated=config["output"] + "output_interomics/vep/cnvs/{sample}_cnvs_vep.vcf",
        sv_annotated=config["output"] + "output_interomics/vep/svs/{sample}_svs_vep.vcf",
    params:
        cache_dir="",  # Diretório do cache (definido no config.yaml)
        species="homo_sapiens"  # Espécie
    conda:
        "../envs/vep.yaml"
    threads: 10
    shell:
        """
        #sed 's/^chr//' {input.cnv} > {input.cnv}.temp ## REMOVER O PREFIXO Chr
        vep -i {input.cnv} -o {output.cnv_annotated} --format vcf --database
        #rm {input.cnv}.temp
        vep -i {input.svs} -o {output.sv_annotated} --format vcf --database
        """
