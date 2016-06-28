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

## Create set of genomes for Fig. 2.
string_genomes = scan("STRING_genomes.txt", what = "character", sep = "\n" )
labs = unlist(strsplit(string_genomes[1], "\t"))
string_genomes = matrix(unlist(strsplit(string_genomes[-1], "\t")), ncol = 11, byrow = T)

string_genomes_c = scan("species.v10.txt", what = "character", sep = "\n" )
labs_c = unlist(strsplit(string_genomes_c[1], "\t"))
string_genomes_c = matrix(unlist(strsplit(string_genomes_c[-1], "\t")), ncol = 4, byrow = T)

string_genomes=cbind(string_genomes, 
			string_genomes_c[match(string_genomes[, 1], string_genomes_c[, 1]), 2])
colnames(string_genomes) = c(labs, labs_c[2])


target_species = c(511145, 212717, 64091)

ds_genomes = string_genomes[(string_genomes[, 12] == "core") |
			(string_genomes[, 9] %in% string_genomes[string_genomes[, 1] %in% target_species, 9]), ]

library(ape)
rRNA = read.dna("STRING_16S_tid.fa", format = "fasta", skip = 0, 
 nlines = 0, comment.char = "#", 
 as.character = T, as.matrix = F)
rRNA_names = names(rRNA)
names(rRNA) = rRNA_names
ds_genomes = ds_genomes[ds_genomes[, 1] %in% names(rRNA), ]

write.table(ds_genomes, "genomes_fig2.txt", sep = "\t", quote = FALSE, row.names = FALSE)



rRNA = lapply(rRNA, paste, collapse = "")
rRNA = unlist(rRNA)
rRNA = rRNA[match(ds_genomes[, 1], names(rRNA) )]

rRNA = paste(">", names(rRNA), "\n", rRNA, sep = "")
write(rRNA, "STRING_16S_fig2.fa")



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
ds_genomes = read.table("genomes_fig2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 
target_species = c(511145, 212717, 64091)
ds_genomes = ds_genomes[!ds_genomes$TID %in% target_species, ]


## Read gene family per genome information into R. These files are obtained by running 
## perl util_get_species_mappings.pl species.mappings.v10.txt
infiles = as.list(paste("families/families_", ds_genomes$TID, sep=""))
families_all = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_all = lapply(families_all, function(x) unlist(strsplit(x[-1], "\t")))
families_all = lapply(families_all, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])
families_all[[1]][1:5]
## [1] "COG2801" "COG0745" "COG0840" "COG1028" "COG0625"

infiles = as.list(paste("families/families_", target_species, sep=""))
families_targets = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_targets = lapply(families_targets, function(x) unlist(strsplit(x[-1], "\t")))
families_targets = lapply(families_targets, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])



## Compute microbe-microbe functional associations, this might take a while.
functional_associations=matrix(ncol = length(families_all), nrow = length(families_targets))
for(i in 1:length(families_targets))
{
	cat(i, "\n")
	functional_associations[i, ] = unlist(lapply(families_all, functional_association, families2 = families_targets[[i]], network = reference_network))
}
functional_associations[, 1:5]
colnames(functional_associations) = ds_genomes$TID
write.table(functional_associations, "associations_fig2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)






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
ds_genomes = read.table("genomes_fig2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 
target_species = c(511145, 212717, 64091)
ds_genomes = ds_genomes[!ds_genomes$TID %in% target_species, ]


## Read gene family per genome information into R. These files are obtained by running 
## perl util_get_species_mappings.pl species.mappings.v10.txt
infiles = as.list(paste("families/families_", ds_genomes$TID, sep=""))
families_all = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_all = lapply(families_all, function(x) unlist(strsplit(x[-1], "\t")))
families_all = lapply(families_all, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])
families_all[[1]][1:5]
## [1] "COG2801" "COG0745" "COG0840" "COG1028" "COG0625"

infiles = as.list(paste("families/families_", target_species, sep=""))
families_targets = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_targets = lapply(families_targets, function(x) unlist(strsplit(x[-1], "\t")))
families_targets = lapply(families_targets, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])



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
sets_all[[1]][[978]]


sets_targets = lapply(families_targets, set_representation, gene_sets = reference_gene_sets)



## Compute genome content similarities, might take a while.
genome_content_similarities=matrix(ncol = length(families_all), nrow = length(families_targets))
for(i in 1:length(families_targets))
{
	cat(i, "\n")
	genome_content_similarities[i, ] = unlist(lapply(sets_all, similarity, gene_sets2 = sets_targets[[i]], threshold = 0.05, size = 4))
}
## Result is simetrical matrix with 1 in diagonal which come from comparing genomes to themselfs.
genome_content_similarities[, 1:5]
colnames(genome_content_similarities) = ds_genomes$TID
write.table(genome_content_similarities, "similarities_fig2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)







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
# ~/Downloads/FastTree_short_branch -nt -gtr STRING_16S_fig2_RDP.fa > STRING_16S_fig2_FastTree




