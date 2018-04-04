
mash_exe=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/mash-Linux64-v2.0/mash
metaBins=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/DNA_BINS/maxbin2_bins
just_seqs=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFERENCE_BINS/just_sequences

#mkdir max_bin_in_thrash
# first pass
#while read a b c d e f g; do
#    echo ${c}
#    $mash_exe screen -p 24 submitted_thrash_genomes.msh ${metaBins}/${c} > max_bin_in_thrash/${c}.tab
#    sort -gr max_bin_in_thrash/${c}.tab > max_bin_in_thrash/${c}.srt.tsv
#    rm max_bin_in_thrash/${c}.tab 
#done < maxbin2_sketch.srted_matches.tsv

# second pass

mkdir thrash_in_maxbin
while read a b c d; do
    echo $b
    $mash_exe screen -p 24 maxbin2_sketch.msh ${just_seqs}/${b} > thrash_in_maxbin/${b}.tab
    sort -gr thrash_in_maxbin/${b}.tab > thrash_in_maxbin/${b}.srt.tsv
    rm thrash_in_maxbin/${b}.tab
done < max_bin_in_thrash/master_containment.tsv

