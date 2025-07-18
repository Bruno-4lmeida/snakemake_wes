rule delly_filter:
    input:
        sv_bcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.bcf",
        cnv_bcf=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.bcf",
    output:
        cnv_vcf=config["output"] + "output_interomics/delly/cnvs/{sample}_cnv.vcf.gz",
        sv_vcf=config["output"] + "output_interomics/delly/svs/{sample}_sv.vcf.gz",
    params:
        reference=config["refgenome"]["ref"],
        mappability=config["refgenome"]["map"]
    conda:
        "../envs/bcftools.yaml"
    threads: 5
    shell:
        """
        echo "##############################################"
        echo "--------   Running Filter Delly  ------------"
        echo "##############################################"

        # Convert CNV BCF to VCF
        bcftools view -i 'FILTER="PASS"' "{input.cnv_bcf}" -Oz -o "{output.cnv_vcf}"
        bcftools index "{output.cnv_vcf}"

        # Convert SV BCF to VCF
        bcftools view -i 'FILTER="PASS"' "{input.sv_bcf}" -Oz -o "{output.sv_vcf}"
        bcftools index "{output.sv_vcf}"
        """
