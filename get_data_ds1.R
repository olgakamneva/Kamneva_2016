############################################################
#### Get OTU-genomes map.
############################################################



setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.
a = scan("clustering_rdp/all_samples.clust", what = "character", sep = "\n")
genomes = scan("STRING_genomes.txt", what = "character", sep = "\n" )
genomes = matrix(unlist(strsplit(genomes[-1], "\t")), ncol = 11, byrow = T)
string16 = scan("STRING_16S_tid.fa", what = "character", sep = "\n" )
string16 = string16[grep(">", string16)]
string16 = sub(">", "", string16)
string16 = matrix(unlist(strsplit(string16, " ")), byrow = T, ncol = 1)



sts = grep("Total Clusters:", a)

for(i in 1:length(sts))
{
	start = sts[i]+1
	if(i == length(sts))
	{
		stop = length(a)
	}
	else
	{
		stop = sts[i+1]-2
	}
	a1 = a[start:stop]
	a11 = matrix(unlist(strsplit(a1, "\t")), ncol = 4, byrow = T)
	a11 = a11[, 4]
	a11_l = strsplit(a11, " ")
	a11_ll = lapply(a11_l, length)
	a11_lgl = lapply(a11_l, function(x) length(grep("OTU", x)))
	a11_ll = unlist(a11_ll)
	a11_lgl = unlist(a11_lgl)
	cat(i, ": ", a[sts[i]-1], length(which(a11_lgl>0 & a11_ll>a11_lgl)), "\n")
}
## 1 :  distance cutoff:	0.0 72 
## 2 :  distance cutoff:	0.01 581 
## 3 :  distance cutoff:	0.02 977 
## 4 :  distance cutoff:	0.03 1119 
## 5 :  distance cutoff:	0.04 1125 
## 6 :  distance cutoff:	0.05 1096 
## 7 :  distance cutoff:	0.06 1017 
## 8 :  distance cutoff:	0.07 941 
## 9 :  distance cutoff:	0.08 851 
## 10 :  distance cutoff:	0.09 740 
## 11 :  distance cutoff:	0.1 645 
## 12 :  distance cutoff:	0.11 544 
## 13 :  distance cutoff:	0.12 459 
## 14 :  distance cutoff:	0.13 394 
## 15 :  distance cutoff:	0.14 332 
## 16 :  distance cutoff:	0.15 288 




i = 4 ## 3% identity cut-off
start = sts[i]+1
stop = sts[i+1]-2
a1 = a[start:stop]
a11 = matrix(unlist(strsplit(a1, "\t")), ncol = 4, byrow = T)
a11 = a11[, 4]
a11_l = strsplit(a11, " ")
a11_ll = lapply(a11_l, length)
a11_lgl = lapply(a11_l, function(x) length(grep("OTU", x)))
a11_ll = unlist(a11_ll)
a11_lgl = unlist(a11_lgl)


keepClusters = which(a11_lgl>0 & a11_ll>a11_lgl)
keepClusters_l = list()
cnt = 0
for(cl in keepClusters)
{
	cnt = cnt+1
	keepClusters_l[[cnt]] = list()
	cluster = a11_l[[cl]]
	
	
	keepClusters_l[[cnt]][[2]] = cluster[grep("OTU", cluster)]
	keepClusters_l[[cnt]][[3]] = cluster[-grep("OTU", cluster)]
	keepClusters_l[[cnt]][[1]] = keepClusters_l[[cnt]][[3]][1]
}
keepClusters_v = c()
for(cl in 1:length(keepClusters_l))
{
	l = length(keepClusters_v)
	keepClusters_v = c(keepClusters_v, paste(keepClusters_l[[cl]][[1]], paste(keepClusters_l[[cl]][[2]], collapse = ", "), paste(keepClusters_l[[cl]][[3]], collapse = ", "), sep = "\t"))
	l1 = length(keepClusters_v)
	if(l1-l>1)
	{
		cat(cl, length(keepClusters_v), "\n")
	}
}

write(keepClusters_v, "gg_string_map")


############################################################
#### Get co-occurrence between relevant OTUs.
############################################################

setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

## 97_99_otu_map.txt obtained by running util_pars_maps.pl which processes 
## 99_otu_map.txt and 99_otu_map.txt file from greengenes to produce 97_99_otu_map.txt.
## 97_99_otu_map.txt contains 97% OTU identifier in column 1 and sequence identifier 
## from 99% OTUs which map onto that 97% OTU in the second column.
OTUmap_v = scan("97_99_otu_map.txt", sep = "\n", what = "character") 
OTUmap_l = strsplit(OTUmap_v, "\t")
OTUmap_t = unlist(OTUmap_l)
OTUmap = matrix(OTUmap_t, ncol = 2, byrow = T)
dim(OTUmap)[1] ## This many sequences are in greengenes
## [1] 1262986


## File gg_13_5_arb_summary_event_OTU was obtained by first running util_pars_arb.pl 
## in the directory where gg_13_5_arb_1.txt, gg_13_5_arb_2.txt, ..., gg_13_5_arb_126.txt, 
## are located to obtain gg_13_5_arb_summary file and then run util_make_event_out_file.pl 
## in the directory where gg_13_5_arb_summary fileis located.
## The file gg_13_5_arb_summary_event_OTU contains gg_id sequence identifier, isolation_source, authors, title fields separated by tabs.

## Files gg_13_5_arb_1.txt, gg_13_5_arb_2.txt, ..., gg_13_5_arb_126.txt and file gg_13_5_arb_summary
## are not provided due to their size.

events_v = scan("gg_13_5_arb_summary_event_OTU", sep = "\n", what = "character")
events_l = strsplit(events_v, "\t")
events_t = unlist(events_l)
events = matrix(events_t, ncol = 4, byrow = T)


head(events)
events = events[-1, ]
events[, 1] = OTUmap[match(events[, 1], OTUmap[, 2]), 1]
events = cbind(events, paste(events[, 2], events[, 3], events[, 4], sep = "_"))
events1 = events



events = events1
head(events)
sps = read.table("gg_string_map", sep = "\t", head = F, colClasses = c("character", "character", "character"))
sps = sps[, c(1, 2)]
sps_l = strsplit(sps[, 2], ", ")
sps = cbind(sps[, 1], unlist(lapply(sps_l, function(x) return(x[1]))))
sps[, 2] = sub("OTU_", "", sps[, 2])
dim(sps)[1] ## This many genomes are mapped onto greengenes.
## [1] 1119
hits = sps
OTUs_my = OTUmap[OTUmap[, 2] %in% hits[, 2], ]


dim(events)[1]  ## Number of greengenes sequences assigned to isolation source.
## [1] 958300

events = events[events[, 1] %in% OTUs_my[, 1], ]
dim(events)[1] ## Number of greengenes sequences assigned to isolation source and mapped onto STRING genomes.
## [1] 53779

events = cbind(events, paste(events[, 1], events[, 2], events[, 3], events[, 4], sep = "_"))
events_l = split(events, as.factor(events[, 6]))
events_l1 = lapply(events_l, matrix, ncol = 6)
events_l2 = lapply(events_l1, function(x) x[1, ])
events2 = matrix(unlist(events_l2), ncol = 6, byrow = T) ## 1 line per OTU per sample
length(unique(ttab[,1]))  ## Number of OTUs
## [1] 1119
length(unique(ttab[,5]))  ## Number of samples
## [1] 6081




events = events2
ttab = events2
s = dim(ttab)[1]
s1 = 0

c = 0;
## Filter to keep OTUs present in at least 2 samples and samples with at least 2 OTUs 
while(s!= s1 & c<20)
{
	c = c+1
	s = dim(ttab)[1]
	ttab_l = split(ttab, as.factor(ttab[, 1]))
	ttab_l_len = unlist(lapply(ttab_l, length))/6
	ttab_l = ttab_l[ttab_l_len>2]
	#ttab_l = ttab_l[ttab_l_len<= 2]
	#ttab_l = ttab_l[1:2]
	ttab_l1 = ttab_l
	ttab_l1_m = lapply(ttab_l1, matrix, ncol = 6)
	ttab_l1_m = lapply(ttab_l1_m, t)
	ttab_l1_m = matrix(unlist(ttab_l1_m), ncol = 6, byrow = T)
	
	ttab = ttab_l1_m
	cat(dim(ttab)[1], "\n")
	ttab_l = split(ttab, as.factor(ttab[, 5]))
	ttab_l_len = unlist(lapply(ttab_l, length))/6
	ttab_l = ttab_l[ttab_l_len>2]
	#ttab_l = ttab_l[ttab_l_len<= 2]
	#ttab_l = ttab_l[1:2]
	ttab_l1 = ttab_l
	ttab_l1_m = lapply(ttab_l1, matrix, ncol = 6)
	ttab_l1_m = lapply(ttab_l1_m, t)
	ttab_l1_m = matrix(unlist(ttab_l1_m), ncol = 6, byrow = T)
	ttab = ttab_l1_m
	
	s1 = dim(ttab)[1]
	cat(dim(ttab)[1], "\n\n")

}

length(unique(ttab[,1]))  ## Number of OTUs. After filterring.
## [1] 308
length(unique(ttab[,5]))  ## Number of samples. After filterring.
## [1] 532




otus = unique(ttab[, 1])
tot = length(unique(ttab[, 5]))
jc_tab = c()
for(i in 1:length(otus))
{
	cat(i, "\n");
	evs1 = ttab[which(ttab[, 1] == otus[i]), 5]
	jc_v = c()
	for(j in 1:length(otus))
	{
		evs2 = ttab[which(ttab[, 1] == otus[j]), 5]	
		jc = length(evs2[evs2 %in% evs1])/length(unique(c(evs1, evs2)))
		jc_v = c(jc_v, jc)
		
	}
	jc_tab = rbind(jc_tab, jc_v)
}

aliase1 = OTUs_my[match(otus, OTUs_my[, 1] ), ]
aliase2 = hits[match(aliase1[, 2], hits[, 2]), ]


string_in_gg = aliase2[, 1]



colnames(jc_tab) = string_in_gg
rownames(jc_tab) = string_in_gg
write.table(jc_tab, "cooccurrence_ds1.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

string_genomes = scan("STRING_genomes.txt", what = "character", sep = "\n" )
labs = unlist(strsplit(string_genomes[1], "\t"))
string_genomes = matrix(unlist(strsplit(string_genomes[-1], "\t")), ncol = 11, byrow = T)
ds_genomes = string_genomes[match(string_in_gg, string_genomes[,1]),]
colnames(ds_genomes)=labs
write.table(ds_genomes, "genomes_ds1.txt", sep = "\t", quote = FALSE, row.names = FALSE)



################################################################################################
#### Compute genome content similarity and microbe-microbe functional 
#### association indices for genomes in dataset 1 using genomics2ecology R package.
################################################################################################

if(!require(devtools, quietly = TRUE))
{
	install.packages(c("devtools"))
}
if(!require(genomics2ecology, quietly = TRUE))
{
	devtools::install_github("olgakamneva/genomics2ecology")
}
library(genomics2ecology)

setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

## Read table of genomes into R.
ds_genomes = read.table("genomes_ds1.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 


## Read gene family per genome information into R. These files are obtained by running 
## perl util_get_species_mappings.pl species.mappings.v10.txt
infiles = as.list(paste("families/families_", ds_genomes$TID, sep=""))
families_all = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_all = lapply(families_all, function(x) unlist(strsplit(x[-1], "\t")))
families_all = lapply(families_all, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])
families_all[[1]][1:5]
## [1] "COG2801" "COG0745" "COG0840" "COG1028" "COG0625"


## Compute microbe-microbe functional associations, this might take a while.
functional_associations=matrix(ncol = length(families_all), nrow = length(families_all))
for(i in 1:length(families_all))
{
	cat(i, "\n")
	functional_associations[i, ] = unlist(lapply(families_all, functional_association, families2 = families_all[[i]], network = reference_network))
}

## Result is simetrical matrix with "NaN" in diagonal which come from division by zero.
functional_associations[1:5, 1:5]
rownames(functional_associations) = ds_genomes$TID
write.table(functional_associations, "associations_ds1.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)






if(!require(devtools, quietly = TRUE))
{
	install.packages(c("devtools"))
}
if(!require(genomics2ecology, quietly = TRUE))
{
	devtools::install_github("olgakamneva/genomics2ecology")
}
library(genomics2ecology)


setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

## Read table of genomes into R.
ds_genomes = read.table("genomes_ds1.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 


## Read gene family per genome information into R. These files are obtained by running 
## perl util_get_species_mappings.pl species.mappings.v10.txt
infiles = as.list(paste("families/families_", ds_genomes$TID, sep=""))
families_all = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_all = lapply(families_all, function(x) unlist(strsplit(x[-1], "\t")))
families_all = lapply(families_all, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])


## genomics2ecology package comes with precomputed sets of 
## functionally linked gene families (gene sets) which are
## storred in reference_gene_sets variable.
## Look at a few the reference gene sets.
data(reference_gene_sets)
reference_gene_sets[[977]]
## [1] "NOG81082"  "NOG82665"  "NOG99413"  "NOG86473"  "NOG168693" "NOG08189" 
## [7] "NOG131456" "NOG160268" "NOG170412"
reference_gene_sets[[978]]
## [1] "COG2801"   "COG3344"   "COG2963"   "COG3547"   "COG2826"   "COG3464"  
## [7] "COG3328"   "NOG102101" "NOG64207" 

## Convert family repertoires for each of the genomes into gene set representation summary
## for each of the gene sets.
sets_all = lapply(families_all, set_representation, gene_sets = reference_gene_sets)
sets_all[[1]][[977]]
## [1] 0 0 0 0 0 0 0 0 0
sets_all[[1]][[978]]
## [1] 1 0 1 0 0 1 1 0 0

## Compute genome content similarities, might take a while.
genome_content_similarities=matrix(ncol = length(families_all), nrow = length(families_all))
for(i in 1:length(families_all))
{
	cat(i, "\n")
	genome_content_similarities[i, ] = unlist(lapply(sets_all, similarity, gene_sets2 = sets_all[[i]], threshold = 0.05, size = 4))
}
## Result is simetrical matrix with 1 in diagonal which come from comparing genomes to themselfs.
genome_content_similarities[1:5, 1:5]
rownames(genome_content_similarities) = ds_genomes$TID
write.table(genome_content_similarities, "similarities_ds1.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)







################################################################################################
#### Get 16S rRNA sequences for relevant genomes.
################################################################################################

setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

## Read table of genomes into R.
ds_genomes = read.table("genomes_ds1.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 

library(ape)
rRNA = read.dna("STRING_16S_tid.fa", format = "fasta", skip = 0, 
 nlines = 0, comment.char = "#", 
 as.character = T, as.matrix = F)
rRNA_names = names(rRNA)
names(rRNA) = rRNA_names

rRNA = lapply(rRNA, paste, collapse = "")
rRNA = unlist(rRNA)
rRNA = rRNA[match(ds_genomes$TID, names(rRNA) )]

rRNA = paste(">", names(rRNA), "\n", rRNA, sep = "")
write(rRNA, "STRING_16S_ds1.fa")
# RDP to get alignment
# ~/Downloads/FastTree_short_branch -nt -gtr STRING_16S_ds1_RDP.fa > STRING_16S_ds1_FastTree






