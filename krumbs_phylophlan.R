if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ggtree")
#BiocManager::install("ggtreeExtra")
#devtools::install_github("YuLab-SMU/ggtree")

{
library(treeio)
#library(devtools)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(ggfittext)
library(ggplot2)
library(tidyverse)
library(ggstance)
#library(gplots)
library(mudata2)
library(ggnewscale)
}

##### ALL PHYLOPHLAN

all_bins_taxa_complete = data.frame(read.csv("data/bin_taxa.tab", sep = "\t"))
all_bins_tax = all_bins_taxa_complete[,c("bin","taxa")]

{
tree = read.newick("data/RAxML_bestTree.all_proteins_refined.tre")
tree = ape::root(tree, node = 209)

#pdf(paste(taxa,"raxml_dendro_nodes_tip.pdf", sep = "_"), width = 12, height = 11)
# GENERATE TREE VIEW WITH NODE NUMBERS, USEFUL FOR ASSIGNING LABELS
pd = ggtree(tree, layout="circular", linetype=1, right=TRUE) + 
  
  # ADD NODE NUMBER LABELS
  geom_text2(aes(label=node), size=1, hjust=-.3) +
  
  # INCREASE EXPANSION INSIDE OF TREE
  scale_x_continuous(expand = c(.05, .05)) + geom_tiplab()

pd = ggtree::rotate(ggtree::rotate(ggtree::rotate(pd, 209), 208), 207)
pd = pd %<+% all_bins_tax + 
  geom_tippoint(aes(color = factor(taxa))) + 
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 10))
#ggtree::rotate(pd, 190)
pd
#dev.off()
}

{
color_alpha = "#004949"
color_gamma = "#490092"
color_bacteroid = "#006DDB"
color_desulfo = "#22CF22"
color_pvc = "#b66dff"
color_spiro = "#ffdf4d"
color_fuso = "#009999"
color_firm = "#920000"
color_fibro = "#252525"
color_weird = "#db6d00"


annotation_pairs = list(
  c(355, "Alphaproteobacteria", color_alpha),
  c(213, "Bacillota", color_firm),
  c(281, "Bacteroidota", color_bacteroid),
  c(356, "Desulfovibrionales", color_desulfo),
  c(204, "Gammaproteobacteria", color_gamma),
  c(153, "Fibrobacterota", color_fibro),
  c(212, "Fusobacteriota", color_fuso),
  c(210, "Spirochaetota", color_spiro),
  c(152, "WOR-3", color_weird),
  c(350, "Verrucomicrobiota", color_pvc)
)

# CREATE TABLE FOR LABEL ASSIGNMENT
hilight_table = data.frame(matrix(ncol = 3, nrow = 0))
colnames(hilight_table) <- c('node_id', 'Taxonomy', 'coloring')

# ADD PAIRS TO TABLE
for (pair in annotation_pairs) {
  hilight_table[nrow(hilight_table) + 1,] = pair
}

# ENSURE NODES ARE NUMBERS RATHER THAN STRINGS
hilight_table$node_id = as.integer(hilight_table$node_id)

# EXPLICITLY ORDER TAXONS (FIXES COLORING BUG)
hilight_table = hilight_table[order(hilight_table$Taxonomy),]
}

{
tree = read.newick("data/RAxML_bestTree.all_proteins_refined.tre")
tree = ape::root(tree, node = 209)
#pdf(paste(taxa,"raxml_dendro.pdf", sep = "_"), width = 12, height = 11)
p = ggtree(tree,layout="fan", open.angle=15, linetype=1, right=TRUE) #+ geom_tiplab(align=TRUE,size=2, color="black")
p2 = ggtree::rotate(ggtree::rotate(p, 209), 208)
p2 = rotate_tree(p2, angle = 90)
p3 = p2 %<+% all_bins_tax + 
        #geom_tippoint(aes(color = factor(taxa))) + 
        #theme(legend.position = "right") + 
        geom_highlight(data = hilight_table, mapping=aes(node=node_id,fill=Taxonomy), extendto = 0.715)+
        scale_size_continuous(range = c(0, 100)) + scale_fill_manual(name = "MAG Taxonomy", values=hilight_table$coloring,  guide = guide_legend(order = 1, ncol = 2))
print(p3)


all_bins_complete = all_bins_taxa_complete[,c("bin", "size", "cgcs", "cgcgenes")]

# 
# my_colors = hilight_table$coloring
# names(my_colors) = hilight_table$Taxonomy
# p4 = p3 + new_scale_fill() + geom_fruit(
#   data=all_bins_complete,
#   geom=geom_col,
#   mapping = aes(x = cgcgenes, y = bin, fill = taxa),
#   axis.params=list(
#     axis = "xy",
#     text.angle = -45,
#     hjust = 0,
#     vjust = 0.5,
#     nbreak = 4, text.size = 1.5
#   ), 
# ) + scale_fill_manual(values=my_colors,
#                         guide = "none")


#print(p4)
}

{
#pdf(paste(taxa,"raxml_qs.pdf", sep = "_"), width = 12, height = 11)
#p4b = p4 + new_scale_fill() + geom_fruit(
#  data=all_bins_complete,
#  geom=geom_col,
#  mapping = aes(x = QualityScore, y= bin, fill=c("green")),
#  axis.params=list(
#    axis = "x",
#    text.angle = -45,
#    hjust = 0,
#    vjust = 0.5,
#    nbreak = 4,
#    text.size = 1.5
#  ), 
#)

#print(p4b)
#dev.off()
}

{
all_bins_cazy = data.frame(read.csv("data/bin_cazy.tab", sep = "\t"))
all_bins_cazy = all_bins_cazy %>% 
  pivot_longer(!bins, names_to = "Enzyme", values_to = "count")
all_bins_cazy$Enzyme = factor(all_bins_cazy$Enzyme, levels=c('Red.Algae.CaZy', 'Brown.Algae.CaZy', 'Green.Algae.CaZy','Sulfatases', 'Peptidoglycanases'))
}


bin_source = data.frame(read.csv("data/bin_source.tab", sep = "\t"))
bin_source$Source = factor(bin_source$Source, levels=c('Fish Gut', 'Bioreactor', 'Seaweed Only Bioreactor'))

taxa="all_taxa"
#svg(paste(taxa,"raxml_seaweed_only_cgcs.svg", sep = "_"), width = 12, height = 11)
p5 = p3 + new_scale_fill() + geom_fruit(
  data=all_bins_cazy,
  geom=geom_bar,
  stat="identity",
  orientation="y",
  mapping = aes(x = count, y = bins, fill = Enzyme),
  axis.params=list(
    axis = "xy",
    title.size = 3,
    vjust = 0,
    hjust = 0,
    text.size = 1.5,
    nbreak = 4,
  ),
) + 
  scale_fill_manual(
    name = "Enzyme Class",
    values=c("indianred1", "lightgoldenrod4", "darkgreen", "gold", "grey"), labels = c("Red Algae CAZy", "Brown Algae CAZy", "Green Algae CAZy","Sulfatases", "Peptidoglycanases"))+
    guides(fill=guide_legend(ncol=2, order=2))

p5
p6 = p5 + new_scale_fill() +  
  geom_fruit(data=bin_source, geom=geom_col,
             mapping=aes(y=bins, x = dummy, fill=Source),pwidth = 0.05
             ) + scale_fill_manual(name = "MAG Source", values = c('Fish Gut' = "#4472C4", 'Bioreactor' = "#A6A6A6", "Seaweed Only Bioreactor" = "pink")) +     guides(fill=guide_legend(ncol=2, order=3))
#+ guides(fill=guide_legend(ncol=3, byrow=TRUE))
svg(paste(taxa,"raxml_no_scfas.svg", sep = "_"), width = 11, height = 7)

p6

dev.off()
#######


svg(paste(taxa,"raxml_scfas.svg", sep = "_"), width = 12, height = 11)

bin_scfa = data.frame(read.csv("data/bin_scfa.tab", sep = "\t"))
bin_scfa = bin_scfa %>% 
  pivot_longer(c(2,3,4), names_to = "SCFA", values_to = "presence")
bin_scfa$presence = factor(bin_scfa$presence, levels=c('Acetate', 'Butyrate', 'Propanoate'))


p7 = p6 + new_scale_fill() +
  geom_fruit(
    data=subset(bin_scfa, !is.na(presence)),
    geom=geom_star,
    mapping=aes(x=SCFA, y=MAG_id, fill=presence, starshape=point),
    size=1.5,
    starstroke=0.01,
    pwidth=0.1,
    inherit.aes = FALSE,
    grid.params = NULL
  ) +
  scale_fill_manual(
    name="Short-chain Fatty Acid Production",
    values=c("red", "darkturquoise", "darkgreen", "white"),
    limits = c('Acetate', 'Butyrate', 'Propanoate'),
    guide=guide_legend(keywidth=1, keyheight=1, order=4,
                       override.aes=list(
                         starshape=c("Acetate"=15,
                                     "Propanoate"=15,
                                     "Butyrate"=15),
                         size=2),             
                       ),
    na.translate=FALSE,
  ) +
  scale_starshape_manual(
    values=c(15),
    guide="none"
  ) +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5.5),
    legend.spacing.y = unit(0.02, "cm")
  )

p7
dev.off()
#########


svg(paste(taxa,"raxml_scfas.svg", sep = "_"), width = 8, height = 8)

bin_sulfrd = data.frame(read.csv("data/bin_sulf_reduct.tab", sep = "\t"))
bin_sulfrd$Sulfate.Reduction  = factor(bin_sulfrd$Sulfate.Reduction, levels=c('Complete Pathway', 'Incomplete Pathway', 'No marker gene'))


p8 = p7 + new_scale_fill() +
  geom_fruit(
    data=subset(bin_sulfrd, !is.na(Sulfate.Reduction)),
    geom=geom_star,
    mapping=aes(y=MAG_id, fill=Sulfate.Reduction),
    starshape = 13,
    na.rm = TRUE,
    size=1.5,
    angle=15,
    starstroke=0.001,
    offset=0.1,
    pwidth=0.1,
    inherit.aes = FALSE,
    grid.params = NULL)+
  scale_fill_manual(
    name="Sulfate Reduction",
    values=c("green", "orange", "white"),
    limits = c('Complete Pathway', 'Incomplete Pathway'),
    guide=guide_legend(keywidth=1, keyheight=1, order=4,
                       override.aes=list(
                         starshape=c('Complete Pathway' = 13, 'Incomplete Pathway' = 13),
                         size=2),             
    ),
    na.translate=FALSE,
  )  +
  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5.5),
    legend.spacing.y = unit(0.02, "cm")
  )
p8
dev.off()

