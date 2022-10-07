rename -v 'R1' 'S1_L001_R1_001' *.fastq.gz
rename -v 'R2' 'S1_L001_R2_001' *.fastq.gz
rename -v 'R3' 'S1_L001_R3_001' *.fastq.gz

# LE-583-KM_ATAC
~/cellranger/cellranger-atac-2.0.0/cellranger-atac count --id=LE-583-KM_ATAC_Analyse_Cellranger \
--reference=~/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--fastqs=. \
--sample=L106799_Track-151941,L107320_Track-151944,L107321_Track-151942,L107322_Track-151943 \
--localcores=8 \
--localmem=64