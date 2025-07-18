rule calculate_proportions:
    input:
        vcf_files=expand(config["output"] + "output_interomics/vcf/filtered/deep_variant/{sample}_final.vcf.gz", sample=SAMPLES),
    output:
        proportions_csv=config["output"] + "output_interomics/metrics/summary/proportions.csv"
    params:
        config["output"] + "output_interomics/metrics/summary/",
    threads: 2
    shell:
        """
        mkdir -p {params}


        echo "Filename;Proportion" > {output.proportions_csv}

        # Itera sobre cada arquivo VCF
        for vcf_file in {input.vcf_files}; do
            if [ -f "$vcf_file" ]; then
                proportion=$(scripts/calculate_proportion.sh "$vcf_file" 10)
                echo "$vcf_file;$proportion" >> {output.proportions_csv}
            fi
        done

        echo "Results saved in {output.proportions_csv}"
        """
