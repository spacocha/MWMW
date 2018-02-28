# MWMW
metagenomic binning, quantification, and analysis with metawrap

TODO:
 - Pull reads that map to 16S gene from libraries, classify, and assign to bins 
   based on combined objective of taxanomic class & abundance vector matches
 - Add Shell scripts from MARCC to scripts

Contents of `Data/` include:
  - `Bin_Taxonomy.xlsx`: Consensus bin taxonomy & scores 
  - `BinningMethodsText.docx`: Description of methods required by GSC 
  - `MAG_Min_Standards_Description.xlsx`: GSC reference document for req'd MAG info
  - `Mystic_MAG_Quality_Stats.xlsx`: Individual MAG stats incl. # unique tRNAs	# rRNA seqs, completeness, contamination, GC,  # contigs/genes/scaffolds, Coding density, GC pct, GC std, Genome size, contig/scaffold length, N50, and more
  - `Bin_Abundances/`
    - `bin_abundance.tab`: matrix of average bin coverage per sample
    - `bin_counts_normed.tsv`: unit normalized within sample & then normalized by library size 
    - `genome_abundance_heatmap.png`: clustered heatmap of `bin_abundance.tab` with assoc. dendograms
    - `sample_read_count.tab`: vector of read counts per sample for normalizing `bin_abundance.tab`
  - `Blob_Plots/`: deeply unhelpful plots of GC% v. coverage of contigs colored in various ways
  - `KEGG_Annotations/:`
    - `AllProteinsRenamed.faa.bz2`: FASTA of ALL protein sequences with bin name & gene id in header (uploaded to KEGG)
    - `Aggregated_Annotations.tsv`: Combined & quality filtered KEGG annotations by GhostKoala & BlastKoala programs
    - `Select_Annotations.tsv`: Same format as `Aggregated_Annotations.tsv` but restricted to KO numbers of interest
    - `Select_Ks_By_Bin.tsv`: a matrix of select KO hit counts by bin
    - `All_HMM_Hits_raw.tsv`: a matrix of db-CAN & metabolic-hmm annotations with gene id & bin info
    - `HMM_counts_by_bin.tsv`: a matrix of hmm hit counts by bin
    - `gene_abundances_raw.tsv`: a matrix of select annotation coverages by sample (sample normalized)
    - `duplicate_clusters.tsv`: output by Salmon as duplicate sequences while calculating coverage (*by_bin.tsv files need correction)
  - `Krona_Plots/`:
    - `mysticLibs_kronagram.html`: assumed taxanomic abundances using QC'd reads as input
    - `final_assembly_kronagram.html`: assumed taxanomic abundances using assembled contigs as input
  - `QUAST_CoAssembly_Stats/`
    - `report.html`: contains co-assembly descriptive statistics
  - `16S_Info/`
    - `16s_annotations.gff`: contains info on bin #, contig id, gene id, and locus, & length
    - `RibosomalRNA.RDP_classes.txt`: contains taxa classes for assembled 16S rRNA
    - `RibosomalRNA.fa`: contains assembled 16S rRNA sequence with bin number of origin & gene id in header
  - `intermediate_files/`: probably unecessary stuff
