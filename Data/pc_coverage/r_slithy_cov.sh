NUM_THREADS=24
b2_idx=/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Batch_Files/MWMW-master/Data/Bin_Fnas_Gffs_Ffns_and_Faas/ref_fastas/Rslithyformis_b2idx
R1=/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/QC_Renamed/H3_PosControlEMAlly_1.fastq
R2=/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/QC_Renamed/H3_PosControlEMAlly_2.fastq

bowtie2 --threads $NUM_THREADS \
-x $b2_idx \
-1 $R1 \
-2 $R2 \
--no-unal \
-S r_slithy.sam

samtools view -F 4 -b r_slithy.sam > r_slithy.bam
samtools sort r_slithy.bam -o r_slithy.sorted.bam
samtools index r_slithy.sorted.bam
bedtools bamtobed -i r_slithy.sorted.bam > r_slithy_reads.bed
python correct_coverage.py r_slithy_reads.bed # keeps pairs if they are aligned within 1kb
python make_genome_bed.py
bedtools coverage -a r_slithy_genome.bed -b r_slithy_reads.pe.bed > r_slithy.cov.bed
