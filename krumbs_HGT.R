if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ggtree")
#BiocManager::install("ggtreeExtra")
#BiocManager::install("ggmsa")
#devtools::install_github("YuLab-SMU/ggtree")
#devtools::install_github("YuLab-SMU/ggmsa")

library(devtools)
library(treeio)

library(ggtree)
library(ggtreeExtra)
library(ggfittext)
library(ggplot2)
library(tidyverse)
library(ggstance)
library(gplots)
library(mudata2)
library(ggrepel)
library(ggnewscale)
library(ggmsa)
library(Biostrings)

all_prot_taxa = data.frame(read.csv("data/all_cazy_prot_HGT.tab", sep = "\t", header = TRUE))
all_prot_taxa$signalp = factor(all_prot_taxa$signalp, c("Y", "N"))

cazy_class = "GH86"
tree = read.newick("data/GH86.tre", node.label = "support")
# SEPARATE BOOTSTRAP VALUES WITH IDENTICAL SEQUENCES (MARKED AS NA BY FASTTREE)
# CHANGE 1 to 1.5 to identify identical sequences separately 
tree@data[["support"]] = replace_na(tree@data[["support"]], 1)
tree = ape::root(tree, node = 112)
alignment = tidy_msa(readAAStringSet("data/GH86_noheader.aln"))
pd = ggtree(tree, layout="rectangular", color = "black", right=TRUE) + 
  
  # ADD NODE NUMBER LABELS
  #geom_text2(aes(label=support), size=2, hjust=-.3) +
  #geom_text2(aes(label=node), size=1, hjust=-.3) +
  
  # INCREASE EXPANSION INSIDE OF TREE
  scale_x_continuous(expand = c(.05, .7))
print(pd)



fill_pvc = "#b66dff"
fill_bacteroid = "#006DDB"
fill_bacill = "#920000"

p2d = ggtree::collapse(ggtree::collapse(ggtree::collapse(ggtree::collapse(ggtree::collapse(ggtree::collapse(pd, node = 80), node = 122), node = 138), node=118), node=113), node=94) + 
  geom_point2(aes(subset=(node==80)), shape=23, size=2, fill=fill_bacteroid) + 
  geom_point2(aes(subset=(node==122)), shape=23, size=2, fill=fill_bacteroid) +
  geom_point2(aes(subset=(node==138)), shape=23, size=2, fill=fill_bacteroid) +
  geom_point2(aes(subset=(node==118)), shape=23, size=2, fill=fill_bacteroid) +
  geom_point2(aes(subset=(node==113)), shape=23, size=2, fill=fill_pvc) + 
  geom_point2(aes(subset=(node==94)), shape=23, size=2, fill=fill_bacteroid)

print(p2d)

#print(p3)
p2 = pd + new_scale_fill() + geom_fruit(
  data = all_prot_taxa,
  geom=geom_col,
  mapping = aes(x = 1, y = prot, fill = source),
  offset = 1.5,
) + scale_fill_manual(name = "Gene Source", values = c('Fish Gut' = "#4472C4", 'Bioreactor' = "#A6A6A6", 'NCBI nr' = 'black'), guide = guide_legend(order = 2))
print(p2)

p2b = p2 +new_scale_fill() + geom_fruit(
  data = all_prot_taxa,
  geom=geom_col,
  mapping = aes(x = 1, y = prot, fill = signalp),
  offset = 0,
) + scale_fill_manual(name = "Signal Peptide", labels = c("Y" = "Exported", "N" = "Not exported"), values = c("#70AD47", "#FFC000"), guide = guide_legend(order = 3))


p3 = p2b %<+% all_prot_taxa +  geom_tiplab(aes(colour = taxonomy), size = 3, align = T, linetype = NULL) + 
  scale_color_manual(name = paste("GH86", "Gene Taxonomy", sep = " "),
                     guide = guide_legend(order = 1),
                     values = c(
                       "Bacillota" = fill_bacill, 
                       "Bacteroidota" = fill_bacteroid,
                       "Verrucomicrobiota" = fill_pvc))

root <- rootnode(tree)
svg("GH86_full_tree.svg")
print(p3)
dev.off()

p5 = p3 + xlim_tree(2) + new_scale_fill() + geom_facet(offset = 3, data = alignment, geom=geom_msa, font = NULL, panel = "msa", border = NA, width = 4)  + guides(fill = "none") + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)
print(p5)


svg("GH86_gene_alignment_collapsed_test.svg")
print(facet_widths(p5, widths = c(3.5, 2.5)))
dev.off()
