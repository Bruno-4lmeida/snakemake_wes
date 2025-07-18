rule test_singularity:
    output:
        "test_output.txt"
    shell:
        """
        singularity exec deepvariant_latest.sif echo "Singularity is working!" > {output}
        """
