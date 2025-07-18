rule annovar:
    input:
        vcf=[
            config["output"] + "output_interomics/vcf/filtered/haplotype_caller/{sample}_final.vcf.gz" if config["use_gatk"] else None,
            config["output"] + "output_interomics/vcf/filtered/strelka/{sample}_final.vcf.gz" if config["use_strelka"] else None,
            config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz" if config["use_deepvariant"] else None,
        ]
    output:
        annotation=[
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annotation_updated.hg38_multianno.csv" if config["use_gatk"] else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annotation_updated.hg38_multianno.csv" if config["use_strelka"] else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annotation_updated.hg38_multianno.csv" if config["use_deepvariant"] else None,
        ]
    threads: config["threads"]["others"]
    params:
        annotation=[
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annotation_updated" if config["use_gatk"] else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annotation_updated" if config["use_strelka"] else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annotation_updated" if config["use_deepvariant"] else None,
        ],
        out=[
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annovar.input" if config["use_gatk"] else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annovar.input" if config["use_strelka"] else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annovar.input" if config["use_deepvariant"] else None,
        ],
        info=[
            config["output"] + "output_interomics/annovar/haplotype_caller/{sample}_annovar.input.INFO" if config["use_gatk"] else None,
            config["output"] + "output_interomics/annovar/strelka/{sample}_annovar.input.INFO" if config["use_strelka"] else None,
            config["output"] + "output_interomics/annovar/deep_variant/{sample}_annovar.input.INFO" if config["use_deepvariant"] else None,
        ],
        db="/home/operator/InterOmics/apps/annovar/humandb",
        build_ver="hg38"
    benchmark:
        config["output"] + "output_interomics/benchmarks/annotated/{sample}_benchmark.txt"
    shell:
        """
        for vcf, info, out, annotation in zip({input.vcf}, {params.info}, {params.out}, {params.annotation}):
            if vcf and info and out and annotation:
                # Conversão do VCF para o formato ANNOVAR
                /home/operator/InterOmics/apps/annovar/convert2annovar.pl -format vcf4 {vcf} -includeinfo -outfile {info}
                /home/operator/InterOmics/apps/annovar/convert2annovar.pl -format vcf4 {vcf} -outfile {out}
                
                # Anotação das variantes usando ANNOVAR
                /home/operator/InterOmics/apps/annovar/table_annovar.pl {out} {params.db} \
                -buildver {params.build_ver} \
                -out {annotation} \
                -remove \
                -protocol refGene,cytoBand,wgRna,1000g2015aug_all,gnomad30_genome,exac03,dbscsnv11,dbnsfp35a,clinvar_20210123,avsnp154 \
                -operation g,r,r,f,f,f,f,f,f,f \
                -nastring . \
                -csvout
        """
