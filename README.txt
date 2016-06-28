## Olga Kamneva. June 27, 2016.
## Not all the files are deposited on GitHub due to their size and availability elsewhere
97_99_otu_map.txt	Mapping of each 97% greengenes OTU to individual sequences
97_otu_map.txt	Mapping of sequences representing 99% greengenes OTUs to 97% OTUs
99_otu_map.txt	Mapping of individual sequences from greengenes 99% OTUs
analysis.R	R code for running analysis (Fig 2, 3, S1 and Table 1)
associations_ds1.txt	Micorbe-microbe functional association indices for organisms in dataset 1
associations_ds2.txt	Micorbe-microbe functional association indices for organisms in dataset 2
associations_fig2.txt	Micorbe-microbe functional association indices for organisms in Fig. 2
clustering_rdp	RDP clustering of 16S rRNA sequences form greengenes and reference genomes from STRING
COG_links_GN_shifted_norm_th275	Network of gene families connected by unit normalized gene neighbor scores and used in MCL clustering 
competition_ds2_full.txt	Metabolic competition indices for gut micro-organisms in Levy and Borenstein 2013
complementarity_ds2_full.txt	Metabolic complementarity indices for gut micro-organisms in Levy and Borenstein 2013
cooccurrence_ds1.txt
cooccurrence_ds2_full.txt	Co-occurrence between gut micro-organisms from Levy and Borenstein 2013
ds1	RDP alignment results of 16S rRNA sequences for organisms in dataset 1
ds2	RDP alignment results of 16S rRNA sequences for organisms in dataset 1
F2A.pdf	Fig. 2 panel A
F2B.pdf	Fig. 2 panel B
families	tables with genome composition for each genome in STRING
fig2	RDP alignment results of 16S rRNA sequences for organisms in Fig 2
Fig3A.pdf	Fig. 3 panel A
Fig3B.pdf	Fig. 3 panel B
Fig3C.pdf	Fig. 3 panel C
Fig3D.pdf	Fig. 3 panel D
FigS1.pdf	Fig. S1 
genomes_ds1.txt List of genomes in dataste 1
genomes_ds2_full_edited.txt	List of gut micro-organisms from Levy and Borenstein 2013 automatically and manually matched to STRING genomes
genomes_ds2_full.txt  List of gut micro-organisms from Levy and Borenstein 2013 automatically matched to STRING genomes
genomes_ds2.txt List of genomes in dataste 2
genomes_fig2.txt List of genomes in Fig. 2
get_data_ds1.R	R code getting data, dataset 1
get_data_ds2.R	R code getting data, dataset 2
get_data_Fig2.R	R code getting data, Fig. 2
gg_genomes_MCL_cluster_coverage	List of genomes from dataset 1 and percent of proteins included into clusters of functionally linked genes
gg_string_map	Mapping of STRING genomes to individual sequences from greengenes representing 97% OTUs
IMG_JGI_gg_genomes	List of genomes from dataset 1 and information about them from IMJ JGI downloaded with NCBI TIDs so some of them are duplicated
Levy_Borenstein_2013_PNAS_sd01.xlsx	Data file from Levy and Borenstein 2013 paper.
mcl_GN275_I40	Clusters of functionally linked genes identified using MCL clustering of the networks of gene familes
README.txt	This file
similarities_ds1.txt	Genome content similarity indices for organisms in dataset 1
similarities_ds2.txt	Genome content similarity indices for organisms in dataset 2
similarities_fig2.txt	Genome content similarity indices for organisms in Fig. 2
STRING_16S_ds1_FastTree	16S rRNA phylogeny, dataset 1
STRING_16S_ds1_RDP.fa	16S rRNA sequence alignment, dataset 1
STRING_16S_ds1.fa	16S rRNA sequences, dataset 1
STRING_16S_ds2_FastTree	16S rRNA phylogeny, dataset 2
STRING_16S_ds2_RDP.fa	16S rRNA sequence alignment, dataset 2
STRING_16S_ds2.fa	16S rRNA sequences, dataset 2
STRING_16S_fig2_FastTree	16S rRNA phylogeny, Fig. 2
STRING_16S_fig2_RDP.fa	16S rRNA sequence alignment, Fig. 2
STRING_16S_fig2.fa	16S rRNA sequences, Fig. 2
STRING_16S_tid.fa	16S rRNA sequences, genomes in STRING database
STRING_genomes.txt	List of genomes in STRING database
util_count_ORFs_by_sp_2.py	Script to prepare gg_protein.aliases.v10.txt and gg_COG.mappings.v10.txt files for util_count_ORFs_by_sp_3.py. Note, location of protein.aliases.v10.txt and COG.mappings.v10.txt which were downloaded from STRING needs to be modified to run
util_count_ORFs_by_sp_3.py	Script to generate gg_genomes_MCL_cluster_coverage file
util_get_species_mappings.pl Script to generate file in families directory, needs species.mappings.v10.txt file from STRING
util_make_event_out_file.pl Generates gg_13_5_arb_summary_event_OTU file from gg_13_5_arb_summary
util_pars_arb.pl	Generates gg_13_5_arb_summary file from gg_13_5_arb_1.txt, ..., gg_13_5_arb_126.txt files from greengenes database
util_pars_maps.pl	Generates 97_99_otu_map.txt file from 97_otu_map.txt and 99_otu_map.txt