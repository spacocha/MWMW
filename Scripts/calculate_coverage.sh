salmon index -k 21 -p 12 -t Annotated_Gene_Seqs.fa -i AGS_Sal_Ind
B_A_S=/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/QC_Renamed
mkdir -p QuantFiles
while read sample R1s R2s; do
    salmon quant --meta -p 12 -i AGS_Sal_Ind --libType IU -1 $B_A_S/$R1s -2 $B_A_S/$R2s -o QuantFiles/${sample}.quant;
done < $B_A_S/samples_names.txt

