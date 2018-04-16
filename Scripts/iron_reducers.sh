#!/bin/bash

#SBATCH
#SBATCH --job-name=blast_iron
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=parallel
#SBATCH --mail-type=END
#SBATCH --mail-user=k.arorawilliams2gmail.com
#SBATCH --error=blst.err
#SBATCH --output=blst.out

module load blast
db_db=../Data/KEGG_Annotations/Iron_Reducers
out_dir=../Data/KEGG_Annotations/Iron_Reducers/match_dir2
mkdir $out_dir

#windowmasker -in $db_db/Rfer_all_proteins.fa -infmt fasta -mk_counts -out $db_db/Rfer_all_proteins.counts
#windowmasker -in $db_db/Rfer_all_proteins.fa -infmt fasta -ustat $db_db/Rfer_all_proteins.counts -outfmt maskinfo_asn1_bin -out $db_db/Rfer_all_proteins.asnb
#makeblastdb -in $db_db/Rfer_all_proteins.fa -dbtype 'nucl' -out $db_db/Rfer_all_proteins -mask_data $db_db/Rfer_all_proteins.asnb

#segmasker -in $db_db/cytochrome_protein.fasta -out $db_db/cytochrome_protein.mask -outfmt maskinfo_asn1_bin 
#makeblastdb -in $db_db/cytochrome_protein.fasta -dbtype 'prot' -out $db_db/cytochrome_protein -mask_data $db_db/cytochrome_protein.mask

ffn_dir=../Data/Bin_Fnas_Gffs_Ffns_and_Faas/bin_untranslated_genes
faa_dir=../Data/Bin_Fnas_Gffs_Ffns_and_Faas/bin_translated_genes

for qq_ff in `ls $ffn_dir/*.ffn`; do
 s_n_ff=`basename $qq_ff`
 echo $s_n_ff;
 tblastx -num_threads 24 -db $db_db/Rfer_all_proteins\
  -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore nident evalue length'\
  -query $qq_ff -out $out_dir/${s_n_ff}.hits -evalue 1E-190 ;
done

for qq_ff in `ls $faa_dir/*.faa`; do
 s_n_ff=`basename $qq_ff`
 echo $s_n_ff
 blastp -query $qq_ff -num_threads 24 -db $db_db/cytochrome_protein\
  -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore nident evalue length'\
  -out $out_dir/${s_n_ff}.hits -evalue 1E-190 ;
done

