# Section: load packages
{
list.of.packages = c('treeio', 'ggtree', 'ggtreeExtra', 'ggstar', 'ggfittext', 'ggplot2',
                 'ggstance', 'ggnewscale', 'mudata2', 'tidyverse')

# Install packages if missing
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
sapply(list.of.packages, library, character.only = TRUE)
}

# Section: data frame preprocessing
{
# Read in node metadata and descriptions
all_bins_taxa_complete = data.frame(read.csv("data/bin_taxa.tab", sep = "\t"))
all_bins_tax = all_bins_taxa_complete[,c("bin","taxa")]

all_bins_complete = all_bins_taxa_complete[,c("bin", "size", "cgcs", "cgcgenes")]
all_bins_cazy = data.frame(read.csv("data/bin_cazy.tab", sep = "\t"))
all_bins_cazy = all_bins_cazy %>% 
  pivot_longer(!bins, names_to = "Enzyme", values_to = "count")
all_bins_cazy$Enzyme = factor(all_bins_cazy$Enzyme, levels=c('Red.Algae.CaZy', 'Brown.Algae.CaZy', 'Green.Algae.CaZy','Sulfatases', 'Peptidogylcanases'))

bin_source = data.frame(read.csv("data/bin_source.tab", sep = "\t"))
bin_source$Source = factor(bin_source$Source, levels=c('Fish Gut', 'Bioreactor'))

bin_scfa = data.frame(read.csv("data/bin_scfa.tab", sep = "\t"))
bin_scfa = bin_scfa %>% 
  pivot_longer(c(2,3,4), names_to = "SCFA", values_to = "presence")
bin_scfa$presence = factor(bin_scfa$presence, levels=c('Acetate', 'Butyrate', 'Propanoate'))
}


# Section: set node colorings
{
# Set colors for phylogenetic clades on tree
  color_alpha = "#006F6C"
  color_gamma = "#8453b5"
  color_bacteroid = "#006DDB"
  color_desulfo = "#22CF22"
  color_pvc = "#d5b2f7"
  color_spiro = "#ffdf4d"
  color_fuso = "#A7EAF7"
  color_firm = "#BC0000"
  color_fibro = "#252525"
  color_weird = "#db6d00"

# Manually place clade starting nodes, with corresponding color
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


# Read in tree and reroot at a specific position  
tree = read.newick("data/RAxML_bestTree.all_proteins_refined.tre")
tree = ape::root(tree, node = 209)

figure = ggtree(tree,layout="fan", open.angle=15, linetype=1, right=TRUE)
figure = ggtree::rotate(ggtree::rotate(figure, 209), 208)
figure = rotate_tree(figure, angle = 90)

figure = figure %<+% all_bins_tax +
        geom_highlight(data = hilight_table, mapping=aes(node=node_id,fill=Taxonomy), extendto = 0.715, alpha = 1, to.bottom = TRUE) +
        scale_size_continuous(range = c(0, 100)) + scale_fill_manual(name = "MAG Taxonomy", values=hilight_table$coloring,  guide = guide_legend(order = 1, ncol = 2)) +
  
  # add enzyme class rings
 new_scale_fill() + geom_fruit(
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
    values=c("indianred1", "lightgoldenrod4", "darkgreen", "gold", "orange"), labels = c("Red Algae CAZy", "Brown Algae CAZy", "Green Algae CAZy","Sulfatases", "Peptidoglycanases"))+
    guides(fill=guide_legend(ncol=2, order=2)) + new_scale_fill() +  
  
  # add MAG source ring
  geom_fruit(data=bin_source, geom=geom_col,
             mapping=aes(y=bins, x = dummy, fill=Source),pwidth = 0.05
             ) + scale_fill_manual(name = "MAG Source", values = c('Fish Gut' = "#4472C4", 'Enrichment' = "#A6A6A6")) +     guides(fill=guide_legend(ncol=2, order=3)) +

  theme(
    legend.background=element_rect(fill=NA),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=5.5),
    legend.spacing.y = unit(0.02, "cm")
  )
figure

svg("figures/phylophlan_figure.svg", width = 8, height = 8)
figure
dev.off()
