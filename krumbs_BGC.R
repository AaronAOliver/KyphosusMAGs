library(ggplot2)
library("stringr") 
library(tidyr)
library(dplyr)
library("RColorBrewer")
library('tidyverse')
library("ggbeeswarm")
library("cowplot")

# Section: load packages
{
  list.of.packages = c( 'ggplot2', 'cowplot', 'ggbeeswarm', 'tidyverse', 'RColorBrewer', 'tidyr', 'dplyr', 'stringr')
  
  # Install packages if missing
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # Load packages
  sapply(list.of.packages, library, character.only = TRUE)
}


bigscape_classes <- c("RiPPs", 
                      "NRPS",
                      "Terpene",
                      "PKSI",
                      "PKSother", 
                      "Saccharides",
                      "PKS-NRP_Hybrids",
                      "Others")                             

my_classes <- c("RiPP", 
                "NRPS",
                "Terpene",
                "PKS", 
                "PKS", 
                "Other",
                "Other",
                "Other")

my_palettes <- c("Greens",
                 "Reds",
                 "Oranges",
                 "Blues",
                 "Greys")

get_random_color <- Vectorize(function(class){
  hue <- my_colors[as.numeric(class)]
  palette  <- brewer.pal(n=9, name = hue)
  sample(palette, 1)})

{
  bgc_df = subset(as.data.frame(read.csv("data/bin_bgc.tab", sep = "\t")), select = c("longBGC", "class", "product")) %>% 
    group_by(class, product) %>%
    tally( )%>%
    mutate(class = factor(class, levels = bigscape_classes), 
           class = my_classes[as.numeric(class)], 
           class = factor(class, levels = unique(my_classes)), 
           color = map_chr(as.numeric(class), 
                           ~sample(brewer.pal(n=9, my_palettes[.x]),1))) %>%
    arrange(class, desc(n)) %>%
    mutate(product  = factor(product, levels = rev(unique(product))), 
           label = str_replace_all(product, ';', ' / '), 
           ymax = cumsum(n), 
           ymin = ymax - n,
           label_position = if_else(n > 20, ymin +(ymax-ymin)/2, as.double(NA_character_)), 
           color = replace(color, product == "terpene", "orange")) 
  
  bgc_df$label[bgc_df$label == "cyclic-lactone-autoinducer"] <- "cyclic-lactone\nautoinducer"
  bgc_df$label[bgc_df$label == "RRE-containing"] <- "RRE-\ncontaining"
  
  barplot <- ggplot(bgc_df, 
                    aes(x = 0,
                        fill = product,
                        label = label)) +
    geom_bar(aes( y = n), 
             position = "stack", stat ="identity") +
    geom_text(aes(y = label_position),
              size = 3.5,
              color = "black") +
    facet_grid( . ~ class, space = "free", switch="x") +
    scale_fill_manual(values = bgc_df$color) +
    scale_x_continuous(breaks = 0) +
    scale_y_continuous(breaks = seq(0, 200, by = 40)) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_blank(), axis.text=element_text(size=12, colour = "black"),
          strip.text.x = element_text(size = 12, colour = "black")) +
    ylab("Number of BGCs") + xlab("")
  barplot
}



svg("figures/bgc_barplot.svg", width = 10, height = 5)
barplot
dev.off()

#####

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
    c("Alphaproteobacteria", color_alpha),
    c("Bacillota", color_firm),
    c("Bacteroidota", color_bacteroid),
    c("Desulfovibrionales", color_desulfo),
    c("Gammaproteobacteria", color_gamma),
    c("Fibrobacterota", color_fibro),
    c("Fusobacteriota", color_fuso),
    c("Spirochaetota", color_spiro),
    c("WOR-3", color_weird),
    c("Verrucomicrobiota", color_pvc)
  )
  
  annotation_pairs = as.data.frame(do.call(rbind,annotation_pairs))
  annotation_pairs$taxonomy = annotation_pairs$V1
  annotation_pairs$color = annotation_pairs$V2
}

{
  bgc_tax_df = subset(as.data.frame(read.csv("data/bin_bgc.tab", sep = "\t")), select = c("longBGC", "class", "taxonomy")) %>% 
    group_by(class, taxonomy) %>%
    tally( )%>%
    mutate(class = factor(class, levels = bigscape_classes), 
           class = my_classes[as.numeric(class)], 
           class = factor(class, levels = unique(my_classes))) %>%
    group_by(class, taxonomy) %>%
    arrange(class, desc(n))
  
  bgc_tax_df <- merge(x=bgc_tax_df,y=annotation_pairs, 
                      by="taxonomy", all.x=TRUE)
  
  barplot_tax <- ggplot(bgc_tax_df,
                        aes(x = 0,
                            fill = taxonomy)) +
    geom_bar(aes( y = n), 
             position = "stack", stat ="identity") +
    facet_grid(. ~ class, space = "free", switch="x") +
    scale_fill_manual(values = unique(bgc_tax_df$color), name = "Source MAG Taxonomy") +
    scale_x_continuous(breaks = 0) + scale_y_reverse(breaks = seq(200, 0, by = -40)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),   strip.background = element_blank(), legend.text=element_text(size = 10),
          strip.text.x = element_blank(), axis.text=element_text(size=12, colour = "black")) +
    ylab("Number of BGCs") + xlab("")
}

svg("figures/bgc_tax_barplot.svg", width = 10, height = 5)
barplot_tax
dev.off()

#############

bgc_violin_df = subset(as.data.frame(read.csv("data/bin_bgc.tab", sep = "\t")), select = c("longBGC", "class", "taxonomy", "completeness","distance")) %>%
  mutate(class = factor(class, levels = bigscape_classes), 
         class = my_classes[as.numeric(class)], 
         class = factor(class, levels = unique(my_classes)))
barplot_violin = ggplot(bgc_violin_df, aes(x=class, y=distance)) +
  geom_quasirandom(aes(x=class, y=distance, fill = completeness), shape = 21, size = 2, color = "transparent", dodge.width = 0.3, varwidth = TRUE) +
  theme_minimal() + 
  geom_hline(yintercept=900, linetype='dotted', col = 'red') +
  ylab("Distance to Nearest GCF") + xlab("") + theme(axis.text=element_text(size=12, colour = "black"))
barplot_violin$labels$fill <- "BGC Quality"

svg("figures/bgc_violin_barplot.svg", width = 7, height = 5)
barplot_violin
dev.off()


{
  g1 = barplot
  g2 = barplot_tax
  g3 = barplot_violin
}

svg("figures/bgc_combined_barplot_bigtext.svg", width = 11, height = 7)
theme_set(theme_minimal())
plot_grid(plot_grid(
  g1 + theme(legend.position = "none"), ggplot(), g3 + theme(legend.position = "none"), get_legend(g3),
  ncol = 4, labels = c("A", "", "B", ""), align = "hv", axis = "b", rel_widths = c(5, 0.5, 3.5, 1)),
  plot_grid(g2 + theme(legend.position = "none"), get_legend(g2), ggplot(), ggplot(),
            ncol = 4, align = "hv", axis = "b", rel_widths = c(5, 2, 1.5, 1.5)), nrow = 2) 
dev.off()

