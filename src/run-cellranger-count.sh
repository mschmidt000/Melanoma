rename -v 'R1' 'S1_L001_R1_001' *.fastq.gz
rename -v 'R2' 'S1_L001_R2_001' *.fastq.gz
rename -v 'R3' 'S1_L001_R3_001' *.fastq.gz

~/cellranger/cellranger-4.0.0/cellranger count --id=LE-583-KM_Analyse_Cellranger \
--transcriptome=~/cellranger/refdata-gex-GRCh38-2020-A \
--fastqs=. \
--sample=L106791_Track-149246,L106792_Track-149249,L106793_Track-149248,L106794_Track-149247 \
--localcores=8 \
--localmem=64
