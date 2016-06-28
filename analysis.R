### Fig S1
setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.
pdf("FigS1.pdf", h=9, w=7)
a=scan("IMG_JGI_gg_genomes", what="character", sep="\n")
a=strsplit(a, "\t")
a=matrix(unlist(a), ncol=213, byrow=T)
colnames(a)=a[1,]
a=a[-1,]
a=a[!duplicated(a[,14]),]

par(mfrow=c(2,1), family="mono")
hist(as.numeric(a[, 167]),xlim=c(0, 100), ylim=c(0,155), breaks=seq(10,65,5), main="", xlab="ORFs in KEGG, %", col="white", border="white")
abline(v=seq(0, 100,20), lty=2, col="grey")
abline(h=seq(0, 200,50), lty=2, col="grey")
hist(as.numeric(a[, 167]), add=T, breaks=seq(10,65,5),col=rgb(0,0,1,0.5), border=rgb(0,0,1,1))
box()
a1=a
mtext("A", side = 3, line = 1, outer = FALSE, at = -10,cex = 2.5)
a=scan("gg_genomes_MCL_cluster_coverage", what="character", sep="\n")
a=strsplit(a, "\t")
a=matrix(unlist(a), ncol=7, byrow=T)
colnames(a)=a[1,]
a=a[-1,]
a=a[match(a1[,14], a[,1]),]
hist(as.numeric(a[,3])/as.numeric(a[,2])*100, xlim=c(0,100), ylim=c(0,155), main="", xlab="ORFs in predicted pathways, %", breaks=seq(30,70,5), col="white", border="white")
abline(v=seq(0, 100,20), lty=2, col="grey")
abline(h=seq(0, 200,50), lty=2, col="grey")
hist(as.numeric(a[,3])/as.numeric(a[,2])*100, col=rgb(255/255,255/255,102/255,0.7),border=rgb(0,1,0,1),add=T, breaks=seq(30,70,5))
hist(as.numeric(a[,5])/as.numeric(a[,2])*100, col=rgb(0,1,0,0.3),border=rgb(0,1,0,1),add=T, breaks=seq(40,85,5))
hist(as.numeric(a[,7])/as.numeric(a[,2])*100, col=rgb(0,153/255,76/255,0.5),border=rgb(0,153/255,76/255,1), add=T, breaks=seq(50,100,5))
legend(-5, 140, c("6 gene families","4 gene families", "2 gene families"), fill=c(rgb(255/255,255/255,102/255,0.7),rgb(0,1,0,0.3),rgb(0,153/255,76/255,0.5) ), box.col="white",border =c(rgb(0,1,0,1),rgb(0,1,0,1),rgb(0,153/255,76/255,1) ), cex=1 );
text("Predicted pathways size of at least:", x=0, y=145, adj=c(0,1))
mtext("B", side = 3, line = 1, outer = FALSE, at = -10,cex = 2.5)
box()
dev.off()




### Fig 2
library(ape)
library(phangorn)

setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.
target_species = c(511145, 212717, 64091)
similarities = read.table("similarities_fig2.txt", sep="\t", head=F)
associations = read.table("associations_fig2.txt", sep="\t", head=F)
sp_tree = read.tree("STRING_16S_fig2_FastTree")
sp_tree = midpoint(sp_tree)
dist = cophenetic.phylo(sp_tree)
x1 = dist[!rownames(dist) %in% target_species, colnames(dist) == target_species[1]]
x2 = dist[!rownames(dist) %in% target_species, colnames(dist) == target_species[2]]
x3 = dist[!rownames(dist) %in% target_species, colnames(dist) == target_species[3]]




y1=as.numeric(similarities[2, match(names(x1), as.character(similarities[1, ]))])
y2=as.numeric(similarities[3, match(names(x2), as.character(similarities[1, ]))])
y3=as.numeric(similarities[4, match(names(x3), as.character(similarities[1, ]))])

pdf("F2A.pdf", h=4, w=4)
layout(matrix(c(2,1,4,3), ncol=2), widths=c(2,1), heights=c(1,2))
par(mar=c(4,4,0,0), family="mono", mgp = c(2, 1, 0))

#plot(1,1,xlim=c(0, 2), ylim=c(.55,1), type="n",xlab="Phylogenetic distance",ylab="Genome\ncontent similaity")
plot(1,1,xlim=c(0, 2), ylim=c(.55,1), type="n",xlab="",ylab="Genome\ncontent similaity", xaxt="n")
abline(h=c(0,.6,.7, .8, .9, 1), col="grey", lty=2)
abline(v=c(0,.5,1, 1.5, 2), col="grey", lty=2)
points(x1,y1, pch=19,col=rgb(1, 0, 0, 0.3), cex=1)
points(x2,y2, pch=19,col=rgb(0, 1, 0, 0.3), cex=1)
points(x3,y3, pch=19, col=rgb(0, 0, 1, 0.3), cex=1)
axis(side = 1,at = c(0, .5,1,1.5, 2),labels = rep("",5), tick = TRUE)


h1=hist(x1, breaks=seq(0, max(x1)+.1, 0.1), plot=F)
h2=hist(x2, breaks=seq(0, max(x2)+.1, 0.1), plot=F)
h3=hist(x3, breaks=seq(0, max(x3)+.1, 0.1), plot=F)


par(mar=c(0,4,0,0), xaxt="n", yaxt="n", bty="n")
plot(1,1 , type="n", ylim=c(0,4),xlim=c(0, 2), xlab="",ylab="")
rect(h1$breaks[1:(length(h1$breaks)-1)],rep(0,length(h1$density)),h1$breaks[2:length(h1$breaks)], h1$density, col=rgb(1, 0, 0, 0.3), border=rgb(1, 0, 0, 0.7))
rect(h2$breaks[1:(length(h2$breaks)-1)],rep(0,length(h2$density)),h2$breaks[2:length(h2$breaks)], h2$density, col=rgb(0, 1, 0, 0.3), border=rgb(0, 1, 0, 0.7))
rect(h3$breaks[1:(length(h3$breaks)-1)],rep(0,length(h3$density)),h3$breaks[2:length(h3$breaks)], h3$density, col=rgb(0, 0, 1, 0.3), border=rgb(0, 0, 1, 0.7))



h1=hist(y1, breaks=seq(0, max(y1)+.1, 0.02), plot=F)
h2=hist(y2, breaks=seq(0, max(y2)+.1, 0.02), plot=F)
h3=hist(y3, breaks=seq(0, max(y3)+.1, 0.02), plot=F)


par(mar=c(4,0,0,0), xaxt="n", yaxt="n", bty="n")
plot(1, 1,  type="n", xlim=c(0,20),ylim=c(.55,1), xlab="",ylab="")
rect(rep(0,length(h1$density)),h1$breaks[1:(length(h1$breaks)-1)], h1$density,h1$breaks[2:length(h1$breaks)], col=rgb(1, 0, 0, 0.3), border=rgb(1, 0, 0, 0.7))
rect(rep(0,length(h2$density)),h2$breaks[1:(length(h2$breaks)-1)], h2$density,h2$breaks[2:length(h2$breaks)], col=rgb(0, 1, 0, 0.3), border=rgb(0, 1, 0, 0.7))
rect(rep(0,length(h3$density)),h3$breaks[1:(length(h3$breaks)-1)], h3$density,h3$breaks[2:length(h3$breaks)], col=rgb(0, 0, 1, 0.3), border=rgb(0, 0, 1, 0.7))


par(mar=c(0,0,0,0), xaxt="n", yaxt="n", bty="n")
plot(1,1,type="n", xlim=c(0,2), ylim=c(-1,6))
#text(c("E. coli", "C. tetani", "Halobacterium sp."), x=rep(0,3), y=c(3,2,1), col=c(rgb(1, 0, 0, 1),rgb(0, 1, 0, 1),rgb(0, 0, 1, 1)), adj=c(0,1), font=4, cex=.8)

dev.off()




y1=as.numeric(associations[2, match(names(x1), as.character(associations[1, ]))])
y2=as.numeric(associations[3, match(names(x2), as.character(associations[1, ]))])
y3=as.numeric(associations[4, match(names(x3), as.character(associations[1, ]))])




pdf("F2B.pdf", h=4, w=4)

layout(matrix(c(2,1,4,3), ncol=2), widths=c(2,1), heights=c(1,2))

par(mar=c(4,4,0,0), family="mono", mgp = c(2, 1, 0),xaxt="l", yaxt="l", bty="o")
plot(1,1,xlim=c(0, 2), ylim=c(0,.4), type="n", xlab="Phylogenetic distance",ylab="Functional\nassociation")
abline(h=c(0,.1,.2, .3, .4), col="grey", lty=2)
abline(v=c(0,.5,1, 1.5, 2), col="grey", lty=2)
points(x1,y1, pch=19,col=rgb(1, 0, 0, 0.3), cex=1)
points(x2,y2, pch=19,col=rgb(0, 1, 0, 0.3), cex=1)
points(x3,y3, pch=19, col=rgb(0, 0, 1, 0.3), cex=1)



h1=hist(x1, breaks=seq(0, max(x1)+.1, 0.1), plot=F)
h2=hist(x2, breaks=seq(0, max(x2)+.1, 0.1), plot=F)
h3=hist(x3, breaks=seq(0, max(x3)+.1, 0.1), plot=F)


par(mar=c(0,4,0,0), xaxt="n", yaxt="n", bty="n")
plot(1,1 , type="n", ylim=c(0,4),xlim=c(0, 2), xlab="",ylab="")
#rect(h1$breaks[1:(length(h1$breaks)-1)],rep(0,length(h1$density)),h1$breaks[2:length(h1$breaks)], h1$density, col=rgb(1, 0, 0, 0.3), border=rgb(1, 0, 0, 0.7))
#rect(h2$breaks[1:(length(h2$breaks)-1)],rep(0,length(h2$density)),h2$breaks[2:length(h2$breaks)], h2$density, col=rgb(0, 1, 0, 0.3), border=rgb(0, 1, 0, 0.7))
#rect(h3$breaks[1:(length(h3$breaks)-1)],rep(0,length(h3$density)),h3$breaks[2:length(h3$breaks)], h3$density, col=rgb(0, 0, 1, 0.3), border=rgb(0, 0, 1, 0.7))

h1=hist(y1, breaks=seq(0, max(y1)+.1, 0.02), plot=F)
h2=hist(y2, breaks=seq(0, max(y2)+.1, 0.02), plot=F)
h3=hist(y3, breaks=seq(0, max(y3)+.1, 0.02), plot=F)


par(mar=c(4,0,0,0), xaxt="n", yaxt="n", bty="n")
plot(1, 1,  type="n", xlim=c(0,20),ylim=c(0,.4), xlab="",ylab="")
rect(rep(0,length(h1$density)),h1$breaks[1:(length(h1$breaks)-1)], h1$density,h1$breaks[2:length(h1$breaks)], col=rgb(1, 0, 0, 0.3), border=rgb(1, 0, 0, 0.7))
rect(rep(0,length(h2$density)),h2$breaks[1:(length(h2$breaks)-1)], h2$density,h2$breaks[2:length(h2$breaks)], col=rgb(0, 1, 0, 0.3), border=rgb(0, 1, 0, 0.7))
rect(rep(0,length(h3$density)),h3$breaks[1:(length(h3$breaks)-1)], h3$density,h3$breaks[2:length(h3$breaks)], col=rgb(0, 0, 1, 0.3), border=rgb(0, 0, 1, 0.7))


par(mar=c(0,0,0,0), xaxt="n", yaxt="n", bty="n")
plot(1,1,type="n", xlim=c(0,2), ylim=c(-1,6))
#text(c("E. coli", "C. tetani", "Halobacterium sp."), x=rep(0,3), y=c(3,2,1), col=c(rgb(1, 0, 0, 1),rgb(0, 1, 0, 1),rgb(0, 0, 1, 1)), adj=c(0,1), font=4 ,cex=.8)

dev.off()







### Fig 3
library(ape)
library(phangorn)
library(vegan)
setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

cooccurence1 = read.table("cooccurrence_ds1.txt", sep="\t", head=F)
cooccurence1 = as.matrix(cooccurence1)
rownames(cooccurence1) = cooccurence1[, 1]
cooccurence1 = cooccurence1[, -1]
colnames(cooccurence1)=rownames(cooccurence1)

sp_tree1=read.tree("STRING_16S_ds1_FastTree")
sp_tree1=midpoint(sp_tree1)
dist1=cophenetic.phylo(sp_tree1)
dist1=dist1[match(rownames(cooccurence1), rownames(dist1)), match(colnames(cooccurence1), colnames(dist1))]

content1=read.table("similarities_ds1.txt", sep="\t", head=F)
content1 = as.matrix(content1)
rownames(content1) = content1[, 1]
content1 = content1[, -1]
colnames(content1)=rownames(content1)

links1=read.table("associations_ds1.txt", sep="\t", head=F)
links1 = as.matrix(links1)
rownames(links1) = links1[, 1]
links1 = links1[, -1]
colnames(links1)=rownames(links1)





content2=read.table("similarities_ds2.txt", sep="\t", head=F)
content2 = as.matrix(content2)
rownames(content2) = content2[, 1]
content2 = content2[, -1]
colnames(content2)=rownames(content2)

links2=read.table("associations_ds2.txt", sep="\t", head=F)
links2 = as.matrix(links2)
rownames(links2) = links2[, 1]
links2 = links2[, -1]
colnames(links2)=rownames(links2)

cooccurence2=read.table("cooccurrence_ds2_full.txt", sep="\t", head=F, stringsAsFactors=F)
cooccurence2 = as.matrix(cooccurence2)
rownames(cooccurence2) = cooccurence2[, 1]
cooccurence2 = cooccurence2[, -1]
colnames(cooccurence2)=rownames(cooccurence2)
class(cooccurence2) = "numeric"
ds2_genomes = read.table("genomes_ds2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 
cooccurence2 = cooccurence2[rownames(cooccurence2) %in% ds2_genomes$name_in_ds2_full, colnames(cooccurence2) %in% ds2_genomes$name_in_ds2_full]
colnames(cooccurence2) = ds2_genomes[match(colnames(cooccurence2), ds2_genomes[, 1]), 2]
rownames(cooccurence2) = ds2_genomes[match(rownames(cooccurence2), ds2_genomes[, 1]), 2]

sp_tree2=read.tree("STRING_16S_ds2_FastTree")
sp_tree2=midpoint(sp_tree2)
dist2=cophenetic.phylo(sp_tree2)
cooccurence2 = cooccurence2[rownames(cooccurence2) %in% rownames(dist2), colnames(cooccurence2) %in% rownames(dist2)]
links2 = links2[rownames(links2) %in% rownames(dist2), colnames(links2) %in% rownames(dist2)]
content2 = content2[rownames(content2) %in% rownames(dist2), colnames(content2) %in% rownames(dist2)]




mantel.partial(cooccurence1, content1, dist1 ,permutations=9999)
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence1, ydis = content1, zdis = dist1,      permutations = 9999) 
## 
## Mantel statistic r: 0.207 
##       Significance: 1e-04 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0308 0.0402 0.0495 0.0581 
## Permutation: free
## Number of permutations: 9999

mantel.partial(cooccurence1, links1, dist1 ,permutations=9999)
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence1, ydis = links1, zdis = dist1,      permutations = 9999) 
## 
## Mantel statistic r: 0.2437 
##       Significance: 1e-04 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0335 0.0428 0.0502 0.0620 
## Permutation: free
## Number of permutations: 9999

mantel.partial(cooccurence2, content2, dist2 ,permutations=9999)
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence2, ydis = content2, zdis = dist2,      permutations = 9999) 
## 
## Mantel statistic r: 0.1954 
##       Significance: 1e-04 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0324 0.0417 0.0506 0.0610 
## Permutation: free
## Number of permutations: 9999


mantel.partial(cooccurence2, links2, dist2 ,permutations=9999)
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence2, ydis = links2, zdis = dist2,      permutations = 9999) 
## 
## Mantel statistic r: 0.06768 
##       Significance: 0.0133 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0379 0.0497 0.0600 0.0714 
## Permutation: free
## Number of permutations: 9999





m1_1=lm(content1[upper.tri(content1)]~dist1[upper.tri(dist1)])
m2_1=lm(links1[upper.tri(links1)]~dist1[upper.tri(dist1)])
m3_1=lm(cooccurence1[upper.tri(cooccurence1)]~dist1[upper.tri(dist1)])

m1_2=lm(content2[upper.tri(content2)]~dist2[upper.tri(dist2)])
m2_2=lm(links2[upper.tri(links2)]~dist2[upper.tri(dist2)])
m3_2=lm(cooccurence2[upper.tri(cooccurence2)]~dist2[upper.tri(dist2)])





pdf("Fig3A.pdf", h=4, w=4.5);
par(mar=c(5,5,1,1), family="mono" ,mgp = c(3, .8, 0))
plot( m1_1$residuals,  m3_1$residuals, pch=19,col=rgb(0,0,1, 0.1), xlim=c(-0.15, 0.25),ylim=c(-.1, 0.9),cex=.7, main="", ylab="Dataset 1\nAdjusted co-occurence",xaxt="n", xlab="", type="n")
abline(v=c(-.1, 0,.1,.2), col="grey", lty=2)
abline(h=c(0, .2,.4,.6,.8), col="grey", lty=2)
points(m1_1$residuals,  m3_1$residuals, pch=19,col=rgb(0,0,1, 0.1), cex=.7)
axis(side = 1,at = c(-.1, 0,.1,.2),labels = rep("", 4), tick = TRUE)
text("A", x=-.15, y=.9, adj=c(0,1), cex=1.8)
text(bquote("Correlation: 0.207"^"**"), x=-.15, y=.8, col="red", adj=c(0,1), font=2)
dev.off()

pdf("Fig3B.pdf", h=4, w=4.5);
par(mar=c(5,5,1,1), family="mono" ,mgp = c(3, .8, 0))
plot(m2_1$residuals,  m3_1$residuals, pch=19,col=rgb(0,0,1, 0.1),ylim=c(-.1, 0.9), xlim=c(-0.2, 0.3), cex=.7, main="",xaxt="n",yaxt="n", ylab="", xlab="", type="n")
abline(v=c(-.2,-.1, 0,.1,.2,.4, .6, .8), col="grey", lty=2)
abline(h=c(0, .2,.4,.6,.8), col="grey", lty=2)
points(m2_1$residuals,  m3_1$residuals, pch=19,col=rgb(0,0,1, 0.1), cex=.7);
axis(side = 1,at = c(-0.2, -0.1,  0, .1, .2, .4, .6, .8),labels = rep("", 8), tick = TRUE)
axis(side = 2,at = c(0,.2,.4, .6, .8),labels = rep("", 5), tick = TRUE)
text("B", x=-.2, y=.9, adj=c(0,1), cex=1.8)
text(bquote("Correlation: 0.2437"^"**"), x=-.2, y=.8, col="red", adj=c(0,1), font=2)
dev.off()

pdf("Fig3C.pdf", h=4, w=4.5);
par(mar=c(5,5,1,1), family="mono", mgp = c(3, .8, 0))
plot(m1_2$residuals,  m3_2$residuals, pch=19,col=rgb(0,0,1, 0.1), xlim=c(-0.15, 0.25),ylim=c(-.1, 0.9),cex=.7,xaxt="n", main="", ylab="Dataset 2\nAdjusted co-occurence", xlab="Adjusted\ngenome content similarity", type="n")
abline(v=c(-.1, 0,.1,.2, .6, .8), col="grey", lty=2)
abline(h=c(0, .2,.4,.6,.8), col="grey", lty=2)
points(m1_2$residuals,  m3_2$residuals, pch=19,col=rgb(0,0,1, 0.1), cex=.7)
axis(side = 1,at = c(-.1, 0,.1,.2),labels = c("-0.1","0","0.1","0.2"), tick = TRUE)
text("C", x=-.15, y=.9, adj=c(0,1), cex=1.8)
text(bquote("Correlation: 0.1954"^"**"), x=-.15, y=.8, col="red", adj=c(0,1), font=2)
dev.off()

pdf("Fig3D.pdf", h=4, w=4.5);
par(mar=c(5,5,1,1), family="mono" ,mgp = c(3, .8, 0))
plot(m2_2$residuals,  m3_2$residuals, pch=19,col=rgb(0,0,1, 0.1),ylim=c(-.1, 0.9), xlim=c(-0.2, 0.3), cex=.7, main="",xaxt="n",yaxt="n", ylab="", xlab="Adjusted microbe-microbe\nfunctional association", type="n")
abline(v=c(-.2,-.1, 0,.1,.2,.4, .6, .8), col="grey", lty=2)
abline(h=c(0, .2,.4,.6,.8), col="grey", lty=2)
points(m2_2$residuals,  m3_2$residuals, pch=19,col=rgb(0,0,1, 0.1), cex=.7);
axis(side = 1,at = c(-0.2, -0.1,  0, .1, .2, .4, .6, .8),labels = c("-0.2", "-0.1","0.0", "0.1", "0.2", "0.4", "0.6", "0.8"), tick = TRUE)
axis(side = 2,at = c(0,.2,.4, .6, .8),labels = rep("", 5), tick = TRUE)
text("D", x=-.2, y=.9, adj=c(0,1), cex=1.8)
text(bquote("Correlation: 0.06768"^"*"), x=-.2, y=.8, col="red", adj=c(0,1), font=2)
dev.off()






###### table 1
library(ape)
library(phangorn)
library(phytools)
setwd("/Users/okamneva/Projects/bacteria_PLoSCB") ## Change this to where all the files are.

content2=read.table("similarities_ds2.txt", sep="\t", head=F)
content2 = as.matrix(content2)
rownames(content2) = content2[, 1]
content2 = content2[, -1]
colnames(content2)=rownames(content2)

links2=read.table("associations_ds2.txt", sep="\t", head=F)
links2 = as.matrix(links2)
rownames(links2) = links2[, 1]
links2 = links2[, -1]
colnames(links2)=rownames(links2)

cooccurence2=read.table("cooccurrence_ds2_full.txt", sep="\t", head=F, stringsAsFactors=F)
cooccurence2 = as.matrix(cooccurence2)
rownames(cooccurence2) = cooccurence2[, 1]
cooccurence2 = cooccurence2[, -1]
colnames(cooccurence2)=rownames(cooccurence2)
class(cooccurence2) = "numeric"

competition2=read.table("competition_ds2_full.txt", sep="\t", head=F, stringsAsFactors=F)
competition2 = as.matrix(competition2)
rownames(competition2) = competition2[, 1]
competition2 = competition2[, -1]
colnames(competition2)=rownames(competition2)
class(competition2) = "numeric"

complementarity2=read.table("complementarity_ds2_full.txt", sep="\t", head=F, stringsAsFactors=F)
complementarity2 = as.matrix(complementarity2)
rownames(complementarity2) = complementarity2[, 1]
complementarity2 = complementarity2[, -1]
colnames(complementarity2)=rownames(complementarity2)
class(complementarity2) = "numeric"



ds2_genomes = read.table("genomes_ds2.txt", sep = "\t", head = T, stringsAsFactors = FALSE) 

cooccurence2 = cooccurence2[rownames(cooccurence2) %in% ds2_genomes$name_in_ds2_full, colnames(cooccurence2) %in% ds2_genomes$name_in_ds2_full]
colnames(cooccurence2) = ds2_genomes[match(colnames(cooccurence2), ds2_genomes[, 1]), 2]
rownames(cooccurence2) = ds2_genomes[match(rownames(cooccurence2), ds2_genomes[, 1]), 2]

competition2 = competition2[rownames(competition2) %in% ds2_genomes$name_in_ds2_full, colnames(competition2) %in% ds2_genomes$name_in_ds2_full]
colnames(competition2) = ds2_genomes[match(colnames(competition2), ds2_genomes[, 1]), 2]
rownames(competition2) = ds2_genomes[match(rownames(competition2), ds2_genomes[, 1]), 2]

complementarity2 = complementarity2[rownames(complementarity2) %in% ds2_genomes$name_in_ds2_full, colnames(complementarity2) %in% ds2_genomes$name_in_ds2_full]
colnames(complementarity2) = ds2_genomes[match(colnames(complementarity2), ds2_genomes[, 1]), 2]
rownames(complementarity2) = ds2_genomes[match(rownames(complementarity2), ds2_genomes[, 1]), 2]


a=cbind(competition2[upper.tri(competition2)], t(competition2)[upper.tri(competition2)])
competition2[upper.tri(competition2)]=apply(a, 1, mean)
competition2=t(competition2)
competition2[upper.tri(competition2)]=apply(a, 1, mean)

a=cbind(complementarity2[upper.tri(complementarity2)], t(complementarity2)[upper.tri(complementarity2)])
complementarity2[upper.tri(complementarity2)]=apply(a, 1, mean)
complementarity2=t(complementarity2)
complementarity2[upper.tri(complementarity2)]=apply(a, 1, mean)
rm(a)




sp_tree2=read.tree("STRING_16S_ds2_FastTree")
sp_tree2=midpoint(sp_tree2)
dist2=cophenetic.phylo(sp_tree2)
cooccurence2 = cooccurence2[rownames(cooccurence2) %in% rownames(dist2), colnames(cooccurence2) %in% rownames(dist2)]
links2 = links2[rownames(links2) %in% rownames(dist2), colnames(links2) %in% rownames(dist2)]
content2 = content2[rownames(content2) %in% rownames(dist2), colnames(content2) %in% rownames(dist2)]
competition2 = competition2[rownames(competition2) %in% rownames(dist2), colnames(competition2) %in% rownames(dist2)]
complementarity2 = complementarity2[rownames(complementarity2) %in% rownames(dist2), colnames(complementarity2) %in% rownames(dist2)]






m1_2=lm(content2[upper.tri(content2)]~dist2[upper.tri(dist2)])
m2_2=lm(links2[upper.tri(links2)]~dist2[upper.tri(dist2)])
m3_2=lm(cooccurence2[upper.tri(cooccurence2)]~dist2[upper.tri(dist2)])
m4_2=lm(competition2[upper.tri(competition2)]~dist2[upper.tri(dist2)])
m5_2=lm(complementarity2[upper.tri(complementarity2)]~dist2[upper.tri(dist2)])




m1_2_matrix=matrix(ncol=dim(content2)[2], nrow=dim(content2)[1])
m1_2_matrix[upper.tri(m1_2_matrix)]=m1_2$residuals
m1_2_matrix[lower.tri(m1_2_matrix)]=m1_2$residuals
diag( m1_2_matrix)=0

m2_2_matrix=matrix(ncol=dim(links2)[2], nrow=dim(links2)[1])
m2_2_matrix[upper.tri(m2_2_matrix)]=m2_2$residuals
m2_2_matrix[lower.tri(m2_2_matrix)]=m2_2$residuals
diag(m2_2_matrix)=0

m3_2_matrix=matrix(ncol=dim(cooccurence2)[2], nrow=dim(cooccurence2)[1])
m3_2_matrix[upper.tri(m3_2_matrix)]=m3_2$residuals
m3_2_matrix[lower.tri(m3_2_matrix)]=m3_2$residuals
diag(m3_2_matrix)=0

m4_2_matrix=matrix(ncol=dim(competition2)[2], nrow=dim(competition2)[1])
m4_2_matrix[upper.tri(m4_2_matrix)]=m4_2$residuals
m4_2_matrix[lower.tri(m4_2_matrix)]=m4_2$residuals
diag(m4_2_matrix)=0

m5_2_matrix=matrix(ncol=dim(complementarity2)[2], nrow=dim(complementarity2)[1])
m5_2_matrix[upper.tri(m5_2_matrix)]=m5_2$residuals
m5_2_matrix[lower.tri(m5_2_matrix)]=m5_2$residuals
diag(m5_2_matrix)=0





regr_31=multi.mantel(m3_2_matrix, m1_2_matrix, nperm=9999)
regr_32=multi.mantel(m3_2_matrix, m2_2_matrix, nperm=9999)
regr_34=multi.mantel(m3_2_matrix, m4_2_matrix, nperm=9999)
regr_35=multi.mantel(m3_2_matrix, m5_2_matrix, nperm=9999)




regr_31$r.squared   
regr_31$fstatistic
regr_31$probF
regr_31$coefficients
regr_31$tstatistic
regr_31$probt


regr_32$r.squared
regr_32$fstatistic
regr_32$probF
regr_32$coefficients
regr_32$tstatistic
regr_32$probt



regr_34$r.squared
regr_34$fstatistic
regr_34$probF
regr_34$coefficients
regr_34$tstatistic
regr_34$probt


regr_35$r.squared
regr_35$fstatistic
regr_35$probF
regr_35$coefficients
regr_35$tstatistic
regr_35$probt



## > regr_31$r.squared   
## [1] 0.03816573
## > regr_31$fstatistic
## [1] 307.4418
## > regr_31$probF
## [1] 0.00010001
## > regr_31$coefficients
##   (intercept)            X1 
## -1.221937e-17  4.764222e-01 
## > regr_31$tstatistic
##   (intercept)            X1 
## -1.173391e-14  1.753402e+01 
## > regr_31$probt
## (intercept)          X1 
##  1.00000000  0.00010001 
## > 
## > 
## > regr_32$r.squared
## [1] 0.004580256
## > regr_32$fstatistic
## [1] 35.65111
## > regr_32$probF
## [1] 0.00160016
## > regr_32$coefficients
##   (intercept)            X1 
## -1.520090e-17  9.251232e-02 
## > regr_32$tstatistic
##   (intercept)            X1 
## -1.434863e-14  5.970855e+00 
## > regr_32$probt
## (intercept)          X1 
##  1.00000000  0.00160016 
## > 
## > 
## > 
## > regr_34$r.squared
## [1] 0.02343459
## > regr_34$fstatistic
## [1] 185.9284
## > regr_34$probF
## [1] 0.00010001
## > regr_34$coefficients
##   (intercept)            X1 
## -2.290616e-17  1.162305e-01 
## > regr_34$tstatistic
##   (intercept)            X1 
## -2.182960e-14  1.363556e+01 
## > regr_34$probt
## (intercept)          X1 
##  1.00000000  0.00010001 
## > 
## > 
## > regr_35$r.squared
## [1] 0.0188562
## > regr_35$fstatistic
## [1] 148.9057
## > regr_35$probF
## [1] 0.00010001
## > regr_35$coefficients
##   (intercept)            X1 
## -1.676237e-17 -2.380720e-01 
## > regr_35$tstatistic
##   (intercept)            X1 
## -1.593724e-14 -1.220269e+01 
## > regr_35$probt
## (intercept)          X1 
##  1.00000000  0.00010001 


mantel.partial(cooccurence2, competition2, dist2 ,permutations=9999)
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence2, ydis = competition2, zdis = dist2,      permutations = 9999) 
## 
## Mantel statistic r: 0.1531 
##       Significance: 1e-04 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0334 0.0429 0.0517 0.0611 
## Permutation: free
## Number of permutations: 9999

mantel.partial(cooccurence2, complementarity2, dist2 ,permutations=9999)
## 
## Partial Mantel statistic based on Pearson's product-moment correlation 
## 
## Call:
## mantel.partial(xdis = cooccurence2, ydis = complementarity2,      zdis = dist2, permutations = 9999) 
## 
## Mantel statistic r: -0.1373 
##       Significance: 1 
## 
## Upper quantiles of permutations (null model):
##    90%    95%  97.5%    99% 
## 0.0301 0.0380 0.0454 0.0524 
## Permutation: free
## Number of permutations: 9999


