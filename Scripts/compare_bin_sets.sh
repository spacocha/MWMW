mash_exe=/home-3/karoraw1@jhu.edu/work/sprehei1/Keith_Files/mash-Linux64-v2.0/mash
parse_exe=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/Scripts/parse_mash_dist.py
metaBins=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/BIN_REASSEMBLY/reassembled_bins
genomes=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFERENCE_BINS/just_sequences
dir_base=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFINED_BINS
mbatBins=$dir_base/metabat2_bins
maxBins=$dir_base/maxbin2_bins
conBins=$dir_base/concoct_bins
out_dir=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFERENCE_BINS/mash_sketches

# make reference sketch
cd $dir2
$mash_exe sketch -p 24 -n -o submitted_thrash_genomes *.fna
mv submitted_thrash_genomes.msh $out_dir

# make all bin set sketches
cd $dir1
$mash_exe sketch -p 24 -n -o mWreassembled_sketch *.fa
mv mW_reassembled_sketch.msh $out_dir
cd $dir3
$mash_exe sketch -p 24 -n -o metabat2_sketch *.fa
mv metabat2_sketch.msh $out_dir
cd $dir4 
$mash_exe sketch -p 24 -n -o maxbin2_sketch *.fa
mv maxbin2_sketch.msh $out_dir
cd $dir5
$mash_exe sketch -p 24 -n -o mWconcoct_sketch *.fa
mv mWconcoct_sketch.msh $out_dir


# I need to calculate min distance & total dispersion for each ref and summary for min distances for an algorithm
cd $out_dir
for i in `ls *_sketch.msh`; do
    $mash_exe dist -p 24 submitted_thrash_genomes.msh $i > ${i}.hits
    python $parse_exe ${i}.hits 
done

# aggregate hits. 

while read a b c d e f g; do
    $mash_exe screen -p 24 submitted_thrash_genomes.msh ${dir3}/${c} > ${c}.tab
    sort -gr ${c}.tab > ${c}.srt.tsv
    rm ${c}.tab
done < metabat2_sketch.srted_matches.tsv
