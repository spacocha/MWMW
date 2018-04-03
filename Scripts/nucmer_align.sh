
_exed_=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/mummer-4.0.0beta2
_qbd_=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/BIN_REASSEMBLY/reassembled_bins
_rbd_=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFERENCE_BINS/just_sequences
_matches_=/home-3/karoraw1@jhu.edu/scratch/THRASH_LIBS/REFERENCE_BINS/mash_sketches/mwR_containment/master_containment.tsv

while read _q_ _r_ _qn_ _rn_; do
  #$_exed_/nucmer --mum $_rbd_/$_r_ $_qbd_/$_q_ -t 24 -p ${_qn_}${_rn_}
  j=$(basename $_qn_$_rn_ .delta)
  $_exed_/delta-filter -l 1000 -q ${_qn_}${_rn_}.delta > ${j}_filter.delta
  $_exed_/show-coords -c -l -L 1000 -r -T ${j}_filter.delta | gzip > ${j}_filter_coords.txt.gz
  echo $j
done < $_matches_

#https://taylorreiter.github.io/2017-08-10-Visualizing-NUCmer-Output/
#j=$(basename $infile .delta)
#delta-filter -l 1000 -q ${infile} > ${j}_filter.delta
#show-coords -c -l -L 1000 -r -T ${j}_filter.delta | gzip > ${j}_filter_coords.txt.gz
