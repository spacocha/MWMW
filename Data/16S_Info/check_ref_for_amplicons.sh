e_coli_idx=/Users/login/Documents/MysticLakeBins/MWMW/Data/Bin_Fnas_Gffs_Ffns_and_Faas/ref_fastas/bt2_idxs/EcoliB_bt2idx
m_aeru_idx=/Users/login/Documents/MysticLakeBins/MWMW/Data/Bin_Fnas_Gffs_Ffns_and_Faas/ref_fastas/bt2_idxs/Maeruginosa_bt2idx
r_slithy_idx=/Users/login/Documents/MysticLakeBins/MWMW/Data/Bin_Fnas_Gffs_Ffns_and_Faas/ref_fastas/bt2_idxs/Rslithyformis_bt2idx

bash bowtie_mod.sh $e_coli_idx e_coli_dbotu.sam 3
bash bowtie_mod.sh $r_slithy_idx r_slithy_dbotu.sam 3
bash bowtie_mod.sh $m_aeru_idx m_aeru_dbotu.sam 3
