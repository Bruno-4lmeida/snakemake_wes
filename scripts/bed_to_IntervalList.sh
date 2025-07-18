gatk BedToIntervalList \
    -I /home/operator/InterOmics/apps/gatk_resources/SureSelectXT_Human_All_Exon_V7.bed \
    -O output.interval_list \
    -SD reference.dict

PicardCommandLine BedToIntervalList \
    I=reference/hg38_Twist_ILMN_Exome_2.5_Panel_annotated.BED \
    O=reference/hg38_Twist_ILMN_Exome_2.5_Panel.interval_list \
    SD=reference/GRCh38_character_genome.dict

rsync -Ap path/origem operator@10.23.5.101:path/destino