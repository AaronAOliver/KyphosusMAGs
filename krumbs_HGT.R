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
library(ggnewscale)

all_prot_taxa = data.frame(read.csv("data/all_cazy_prot_HGT.tab", sep = "\t", header = TRUE))
all_prot_taxa$signalp = factor(all_prot_taxa$signalp, c("Y", "N"))

cazy_class = "GH86"
tree = read.newick("data/GH86.tre", node.label = "support")
# SEPARATE BOOTSTRAP VALUES WITH IDENTICAL SEQUENCES (MARKED AS NA BY FASTTREE)
# CHANGE 1 to 1.5 to identify identical sequences separately 
tree@data[["support"]] = replace_na(tree@data[["support"]], 1)
tree = ape::root(tree, node = 94)
alignment = tidy_msa(readAAStringSet("data/GH86_noheader.aln"))
pd = ggtree(tree, layout="rectangular", color = "black", right=TRUE)
  
  # ADD NODE NUMBER LABELS
  #geom_text2(aes(label=support), size=2, hjust=-.3) +
  #geom_text2(aes(label=node), size=2, hjust=-.3)
  
  # INCREASE EXPANSION INSIDE OF TREE
  #scale_x_continuous(expand = c(.1, .7))
print(pd)



fill_pvc = "#b66dff"
fill_bacteroid = "#006DDB"
fill_bacill = "#920000"

p2d = ggtree::collapse(ggtree::collapse(pd, node = 94), node = 112)+ 
  geom_point2(aes(subset=(node==112)), shape=23, size=2, fill="black") + 
  geom_point2(aes(subset=(node==94)), shape=23, size=2, fill="black")

print(p2d)

#print(p3)
p2 = p2d + new_scale_fill() + geom_fruit(
  data = all_prot_taxa,
  geom=geom_col,
  mapping = aes(x = 1, y = prot, fill = source),
  offset = 1.5,
) + scale_fill_manual(name = "Gene Source", values = c('Fish Gut' = "#4472C4", 'Bioreactor' = "#A6A6A6", 'NCBI nr' = 'black'), guide = guide_legend(order = 2))
print(p2)

p2b = p2d +new_scale_fill() + geom_fruit(
  data = all_prot_taxa,
  geom=geom_col,
  mapping = aes(x = 1, y = prot, fill = signalp),
  offset = 0,
) + scale_fill_manual(name = "Signal Peptide", labels = c("Y" = "Exported", "N" = "Not exported"), values = c("#70AD47", "#FFC000"), guide = guide_legend(order = 3))


p3 = p2d %<+% all_prot_taxa +  geom_tiplab(aes(colour = taxonomy), size = 3, align = T, linetype = NULL) + 
  scale_color_manual(name = paste("GH86", "Gene Taxonomy", sep = " "),
                     guide = guide_legend(order = 1),
                     values = c(
                       "Bacillota" = fill_bacill, 
                       "Bacteroidota" = fill_bacteroid,
                       "Verrucomicrobiota" = fill_pvc))







################ Alignment figure
root <- rootnode(tree)
#svg("GH86_full_tree.svg")
print(p3)
#dev.off()

{
get_monocolor_coloring = function(coloring="black"){
  data.frame(names = c(LETTERS[1:26],"-"), 
             color = c(rep(coloring,each=26), "white"), 
             stringsAsFactors = FALSE)
}

signal_cutoff = 45
undescribed_cutoff = 225
catalytic_cutoff = 1060
total = 1296

aln_signal = alignment[alignment$position <= signal_cutoff, ]
aln_undescribed = alignment[alignment$position > signal_cutoff & alignment$position <= undescribed_cutoff, ]
aln_catalytic = alignment[alignment$position > undescribed_cutoff,  ]

mycustom_signal <- get_monocolor_coloring("#70AD47")
mycustom_undescribed <- get_monocolor_coloring("#AA336A")
mycustom_catalytic <- get_monocolor_coloring("#56494e")
}

p5 = p3 + xlim_tree(2) +
  new_scale_fill() + geom_facet(offset = 0, data = aln_signal, geom=geom_msa, font = NULL, panel = "signal", border = NA, width = 0.0925, custom_color = mycustom_signal)  + 
  new_scale_fill() + geom_facet(offset = 0, data = aln_undescribed, geom=geom_msa, font = NULL, panel = "undescribed", border = NA, width = 4, custom_color = mycustom_undescribed) +
  new_scale_fill() + geom_facet(offset = 0, data = aln_catalytic, geom=geom_msa, font = NULL, panel = "catalytic", border = NA, width = 4, custom_color = mycustom_catalytic) +
  guides(fill = "none") + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  panel.spacing.x = unit(0, "pt")
  
) + xlim_tree(6)
#print(p5)


summed_width = 5
svg("figures/GH86_gene_alignment.svg", width = 6, height = 4)
print(facet_widths(p5, widths = c(3, summed_width * signal_cutoff/total,summed_width * (undescribed_cutoff - signal_cutoff)/total, 2 * summed_width * (total - catalytic_cutoff)/total)))
dev.off()


################# For CMBB talk

all_prot_taxa = data.frame(read.csv("data/all_cazy_prot_HGT.tab", sep = "\t", header = TRUE))
all_prot_taxa$signalp = factor(all_prot_taxa$signalp, c("Y", "N"))

cazy_class = "GH86"
tree = read.newick("data/GH86.tre", node.label = "support")
# SEPARATE BOOTSTRAP VALUES WITH IDENTICAL SEQUENCES (MARKED AS NA BY FASTTREE)
# CHANGE 1 to 1.5 to identify identical sequences separately 
tree@data[["support"]] = replace_na(tree@data[["support"]], 1)
tree = ape::root(tree, node = 112)

pd = ggtree(tree, layout="circular", color = "black", right=TRUE) + 
  
  # ADD NODE NUMBER LABELS
  #geom_text2(aes(label=support), size=2, hjust=-.3) +
  #geom_text2(aes(label=node), size=1, hjust=-.3) +
  
  # INCREASE EXPANSION INSIDE OF TREE
  scale_x_continuous(expand = c(.05, .7))
print(pd)




# Section: set node colorings
{
  # Set colors for phylogenetic clades on tree
  color_pvc = "#b66dff"
  color_bacteroid = "#006DDB"
  color_firm = "#920000"
  
  
  # Manually place clade starting nodes, with corresponding color
  annotation_pairs = list(
    c(110, "Bacillota", color_firm),
    c(117, "Bacteroidota", color_bacteroid),
    c(94, "Bacteroidota", color_bacteroid),
    c(98, "Bacteroidota", color_bacteroid),
    c(109, "Verrucomicrobiota", color_pvc),
    c(105, "Verrucomicrobiota", color_pvc),
    c(113, "Verrucomicrobiota", color_pvc)
  )
  
  # Create table for taxonomy assignment
  hilight_table = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(hilight_table) <- c('node_id', 'Taxonomy', 'coloring')
  
  # Add taxonomy pairs to table
  for (pair in annotation_pairs) {
    hilight_table[nrow(hilight_table) + 1,] = pair
  }
  
  # Ensure nodes are numbers rather than strings
  hilight_table$node_id = as.integer(hilight_table$node_id)
  
  # Explicitly order taxonomic groups
  hilight_table = hilight_table[order(hilight_table$Taxonomy),]
}

figure = pd %<+% all_prot_taxa + 
  geom_highlight(data = hilight_table, mapping=aes(node=node_id,fill=Taxonomy),  to.bottom = TRUE)+
  scale_size_continuous(range = c(0, 100)) + scale_fill_manual(name = "Gene Taxonomy", values = c(
    "Bacillota" = color_firm, 
    "Bacteroidota" = color_bacteroid,
    "Verrucomicrobiota" = color_pvc),  guide = guide_legend(order = 1, ncol = 1))

print(figure)

ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

p3 = pd %<+% all_prot_taxa + 
  geom_highlight(data = hilight_table, mapping=aes(node=node_id,fill=Taxonomy), extendto = 0.715) + 
  geom_tiplab(aes(colour = taxonomy), size = 3, align = T, linetype = NULL) + 
  scale_color_manual(name = paste("GH86", "Gene Taxonomy", sep = " "),
                     guide = guide_legend(order = 1),
                     values = c(
                       "Bacillota" = fill_bacill, 
                       "Bacteroidota" = fill_bacteroid,
                       "Verrucomicrobiota" = fill_pvc))

svg("GH86_draft_figure")
print(figure)
dev.off()
