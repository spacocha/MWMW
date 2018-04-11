
module load blast
db_db=../Data/KEGG_Annotations/Iron_Reducers
out_dir=../Data/KEGG_Annotations/Iron_Reducers/match_dir
mkdir $out_dir
makeblastdb -in $db_db/Rfer_all_proteins.fa -dbtype 'nucl' -out $db_db/Rfer_all_proteins 

ffn_dir=../Data/Bin_Fnas_Gffs_Ffns_and_Faas/bin_untranslated_genes
for qq_ff in `ls $ffn_dir`; do
echo $qq_ff;
tblastx -num_threads 24 -db $db_db/Rfer_all_proteins -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore nident evalue length' -query $ffn_dir/$qq_ff -out $out_dir/${qq_ff}.hits -evalue 1E-150;
#Add something here to process all of the results into a table
done

makeblastdb -in $db_db/cytochrome_protein.fasta -dbtype 'prot' -out $db_db/cytochrome_protein 


faa_dir=../Data/Bin_Fnas_Gffs_Ffns_and_Faas/bin_translated_genes
for qq_ff in `ls $faa_dir`; do
echo $qq_ff;
blastp -query $faa_dir/$qq_ff -num_threads 24 -db $db_db/cytochrome_protein -outfmt '6 qseqid qstart qend qlen sseqid staxids sstart send bitscore nident evalue length' -out $out_dir/${qq_ff}.hits -evalue 1E-150;
#Add something here to process all the results into a table
done

