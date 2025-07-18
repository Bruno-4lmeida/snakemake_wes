rule sansa_vcf:
    input:
        cnv_anno=config["output"] + "output_interomics/sansa/cnvs/{sample}_cnv_anno.bcf",
        sv_anno=config["output"] + "output_interomics/sansa/svs/{sample}_SV_anno.bcf",
    output:
        cnv_vcf=config["output"] + "output_interomics/sansa/cnvs/{sample}_cnv_anno.vcf.gz",
        sv_vcf=config["output"] + "output_interomics/sansa/svs/{sample}_SV_anno.vcf.gz",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    shell:
        """
        echo "##############################################"
        echo "---- Converting bcf sansa to VCF.gz ---------"
        echo "##############################################"

        bcftools view "{input.cnv_anno}" -Oz -o "{output.cnv_vcf}"
        bcftools index "{output.cnv_vcf}"

        # Convert SV BCF to VCF
        bcftools view "{input.sv_anno}" -Oz -o "{output.sv_vcf}"
        bcftools index "{output.sv_vcf}"

        """
