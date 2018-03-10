bowtie2 --score-min 'C,0,-1' --no-sq -f --all --end-to-end --threads $3 -x $1 \
-U unique.dbOTU.nonchimera.ns.fasta --no-unal -S $2

