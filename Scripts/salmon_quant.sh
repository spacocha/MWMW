#!/bin/bash

#SBATCH
#SBATCH --job-name=salmon_quant
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G
#SBATCH --partition=gpu
#SBATCH --mail-type=END
#SBATCH --mail-user=k.arorawilliams2gmail.com
#SBATCH --error=sq.err
#SBATCH --output=sq.out

source activate metawrap2-env
sq_dir=/home-3/karoraw1@jhu.edu/scratch/MWMW/Data/KEGG_Annotations
python select_gene_seqs.py $sq_dir
salmon index -k 21 -p 12 -t $sq_dir/Annotated_Gene_Seqs.fa -i $sq_dir/AGS_Sal_Ind

B_A_S=/home-3/karoraw1@jhu.edu/scratch/metaWRAP_Out/QC_Renamed

mkdir -p $sq_dir/QuantFiles

while read sample R1s R2s; do
    salmon quant --meta -i $sq_dir/AGS_Sal_Ind --libType IU -1 $B_A_S/$R1s -2 $B_A_S/$R2s -o $sq_dir/QuantFiles/${sample}.quant;
done < $B_A_S/samples_names.txt

sc_dir=`pwd`
cd $sq_dir/QuantFiles
python summarize_salmon.py
python $sc_dir/concatenate_salmon.py $sq_dir
