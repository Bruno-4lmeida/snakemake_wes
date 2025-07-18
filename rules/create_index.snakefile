configfile: "config.yaml"

rule create_index:
    output:
        fasta=config["refgenome"]["ref_dir"] + "GRCh38_genome.fasta",
        dict=config["refgenome"]["ref_dir"] + "GRCh38_genome.dict",
        fai=config["refgenome"]["ref_dir"] + "GRCh38_genome.fasta.fai",
        bwa_index=config["refgenome"]["ref_dir"] + "GRCh38_genome.fasta.bwt"
    conda:
        "../envs/index.yaml", 
    shell:
        """
        echo "Baixando o genoma de referência..."
        wget -O {output.fasta}.gz ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        
        echo "Descompactando o genoma..."
        gunzip -c {output.fasta}.gz > {output.fasta}

        echo "Criando Sequence Dictionary..."
        PicardCommandLine CreateSequenceDictionary R={output.fasta} O={output.dict}

        echo "Criando arquivo .fai..."
        samtools faidx {output.fasta}

        echo "Criando índice BWA..."
        bwa index {output.fasta}
        """
