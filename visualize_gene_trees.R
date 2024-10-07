
library(ape)

### Reading in gene trees ------
## "all_genes.tre" in the below line is a placeholder
## you should give a name and path of a file that has all your gene trees
## you can create it by "cat"ing all your *.tre files from you genes
gene_trees = read.tree("all_trees_good.tre") 

## I alsa had written the busco IDs/gene names of all the gene trees that are
## present in the above "all_genes.tre" file into another file called
## "gene_list.txt" this file should have the names in the same order as the
## gene trees appear in your "all_genes.tre". This is just to put titles for
## each gene tree plot so that I would know which gene I am looking at.
gene_names = scan("name_final.txt", what=character())
names(gene_trees) = gene_names

### checking the remaining 'ok' trees and hybpiper candidate paralogs
#paralogs = scan("all_paralogs_name.txt", character())
#paralogs = paste(paralogs, ".paralogs.tre", sep="")
##non_weird_paralogs = setdiff(paralogs, names(weird_genes))

## nsamples = sapply(gene_trees, function(x){length(x$tip.label)})
## names(nsamples) = names(gene_trees)
## tmp = nsamples[nsamples > 66]
## length(tmp)

#pdf("paralogs.pdf", 30,30)
#par(mfrow=c(5,5))
#sapply(paralogs, function(x){
  ##print(x)
#  plot(gene_trees[[x]], main=x, cex=0.5)
#})
#dev.off()


## Plotting the relative branch lengths of the longest branches

brlens = lapply(gene_trees, function(x){x$edge.len})
prop_max = sapply(brlens, function(x){max(x)/sum(x)})
names(prop_max) = gene_names

#prop_max = prop_max[setdiff(names(prop_max), paralogs)]

plot(density(prop_max, na.rm = T))
hist(prop_max, br=40)

threshold = 0.03 ## upper limit for the trees that regarded as 'normal'
weird_genes = prop_max[prop_max > threshold]

pdf("good_genes_fsa.pdf", 30,30)
par(mfrow=c(5,5))
sapply(names(weird_genes), function(x){
 # print(x)
  plot(gene_trees[[x]], main=x)
})
dev.off()

barplot(table(sapply(setdiff(gene_names, names(weird_genes)), function(x){
  length(gene_trees[[x]]$tip.label)
})), las=2)

write(names(weird_genes), file = "outlier_genes_fsa_005.txt",sep='\n')

genes_to_remove = union(names(weird_genes), paralogs)


nsamples = sapply(gene_trees, function(x){length(x$tip.label)})

nsamples_ok_genes = nsamples[setdiff(names(nsamples), genes_to_remove)]

par(mfrow=c(1,2))
barplot(table(nsamples), las=2)
barplot(table(nsamples_ok_genes), las=2)

par(mfrow=c(1,1))
prop_samples = nsamples_ok_genes / 66
quantile(prop_samples)
hist(prop_samples)
abline(v=c(0.6,0.7,0.8))

sum(prop_samples >= 0.7)
