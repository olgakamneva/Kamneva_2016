###########################################################################
#### Match genomes from Levy, Borenstein. 2013. PNAS to 
#### STRING genomes to create dataset 2.
###########################################################################


setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

string_genomes = scan("STRING_genomes.txt", what = "character", sep = "\n" )
string_genomes = matrix(unlist(strsplit(string_genomes[-1], "\t")), ncol = 11, byrow = T)
string_labels = string_genomes[, 3]
string_labels = gsub(" |-|/|\\.", "_", string_labels)
string_labels = gsub("__", "_", string_labels)

ds = read.table("competition_ds2_full.txt", sep = "\t", head = F)
dim(ds)[1]
## [1] 154  -- This many genomes are in original file.
ds_labels = as.character(ds[, 1])

length(string_labels[string_labels %in% ds_labels])
## [1] 83  -- This many matched exactly.

ds_genomes = cbind(ds_labels, string_genomes[match(ds_labels, string_labels),])
write.table(ds_genomes, "genomes_ds2_full.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ds_genomes = read.table("genomes_ds2_full_edited.txt", sep = "\t", head = T, stringsAsFactors = FALSE) ## Input file was edited by hand as described in the manuscript.
ds_genomes = ds_genomes[!is.na(ds_genomes[,2]), ]
dim(ds_genomes)[1]
## [1] 128 -- This many were matched in total and are in dataset 2.
write.table(ds_genomes, "genomes_ds2.txt", sep = "\t", quote = FALSE, row.names = FALSE)



###########################################################################
#### Compute genome content similarity and microbe-microbe functional 
#### association indices for genomes using genomics2ecology R package.
###########################################################################

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
ds_genomes = read.table("genomes_ds2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 


## Read gene family per genome information into R. These files are obtained by running 
## perl util_get_species_mappings.pl species.mappings.v10.txt
infiles = as.list(paste("families/families_", ds_genomes$TID, sep=""))
families_all = lapply(infiles, scan, what = "character", sep = "\n", quiet = TRUE)
families_all = lapply(families_all, function(x) unlist(strsplit(x[-1], "\t")))
families_all = lapply(families_all, function(x) matrix(x, ncol = 3, byrow = TRUE)[,2])
families_all[[1]][1:5]
## [1] "COG2801"  "NOG64207" "COG0582"  "NOG00011" "COG3436"  

## Compute microbe-microbe functional associations, this might take a while.
functional_associations=matrix(ncol = length(families_all), nrow = length(families_all))
for(i in 1:length(families_all))
{
	cat(i, "\n")
	functional_associations[i, ] = unlist(lapply(families_all, functional_association, families2 = families_all[[i]], network = reference_network))
}

## Result is simetrical matrix with "NaN" in diagonal which come from division by zero.
functional_associations[1:5, 1:5]
##            [,1]      [,2]      [,3]       [,4]       [,5]
## [1,]        NaN 0.1017190 0.1187258 0.08891993 0.06871782
## [2,] 0.10171896       NaN 0.1971182 0.30640569 0.24045352
## [3,] 0.11872580 0.1971182       NaN 0.17684315 0.13661285
## [4,] 0.08891993 0.3064057 0.1768431        NaN 0.20227671
## [5,] 0.06871782 0.2404535 0.1366129 0.20227671        NaN

rownames(functional_associations) = ds_genomes$TID
write.table(functional_associations, "associations_ds2.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)








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
## [1] 1 1 1 0 1 0 1 0 1

## Compute genome content similarities, might take a while.
genome_content_similarities=matrix(ncol = length(families_all), nrow = length(families_all))
for(i in 1:length(families_all))
{
	cat(i, "\n")
	genome_content_similarities[i, ] = unlist(lapply(sets_all, similarity, gene_sets2 = sets_all[[i]], threshold = 0.05, size = 4))
}
## Result is simetrical matrix with 1 in diagonal which come from comparing genomes to themselfs.
genome_content_similarities[1:5, 1:5]
##           [,1]      [,2]      [,3]      [,4]      [,5]
## [1,] 1.0000000 0.6924515 0.7025868 0.6980224 0.6983930
## [2,] 0.6924515 1.0000000 0.7584959 0.7349576 0.7434696
## [3,] 0.7025868 0.7584959 1.0000000 0.7343558 0.7544916
## [4,] 0.6980224 0.7349576 0.7343558 1.0000000 0.7593053
## [5,] 0.6983930 0.7434696 0.7544916 0.7593053 1.0000000


rownames(genome_content_similarities) = ds_genomes$TID
write.table(genome_content_similarities, "similarities_ds2.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)




################################################################################################
#### Get 16S rRNA sequences for relevant genomes.
################################################################################################

setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

## Read table of genomes into R.
ds_genomes = read.table("genomes_ds2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 

library(ape)
rRNA = read.dna("STRING_16S_tid.fa", format = "fasta", skip = 0, 
 nlines = 0, comment.char = "#", 
 as.character = T, as.matrix = F)
rRNA_names = names(rRNA)
names(rRNA) = rRNA_names
ds_genomes = ds_genomes[ds_genomes$TID %in% names(rRNA), ]

rRNA = lapply(rRNA, paste, collapse = "")
rRNA = unlist(rRNA)
rRNA = rRNA[match(ds_genomes$TID, names(rRNA) )]

rRNA = paste(">", names(rRNA), "\n", rRNA, sep = "")
write(rRNA, "STRING_16S_ds2.fa")
# RDP to get alignment
# ~/Downloads/FastTree_short_branch -nt -gtr STRING_16S_ds2_RDP.fa > STRING_16S_ds2_FastTree






