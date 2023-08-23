# Get visualization packages
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('ggbeeswarm')) install.packages('ggbeeswarm'); library('ggbeeswarm')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('stringr')) install.packages('stringr'); library('stringr')

# Get data manipulation packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')

# Read data from file
data <- read.csv("data/read_taxonomy.tab", stringsAsFactors = FALSE, sep="\t")
sample_order = c("R1_6", "R1_12", "R1_24", "R1_25", "R1_28", "R1_29", "R1_38", "R1_43", "R2_8", "R2_24", "R2_26", "Empty",
                 "F5GI_2", "F5GI_3", "F5HG_2", "F5HG_3", "F6GI_3", "F6GI_4", "F6HG_2", "F6HG_3",
                 "F7GI_2", "F7GI_3", "F7HG_2", "F7HG_3", "F8GI_2", "F8GI_3", "F8HG_2", "F8HG_3")
data$id = factor(data$id, levels = sample_order)

taxa_order = rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                   "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                 "WOR.3", "Other.Taxa"))

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
color_euk = "#696969"
color_other = "#D3D3D3"
}

# Pivot data using tidyverse
data_pivoted <- data %>%
  pivot_longer(
    cols = taxa_order,
    names_to = "taxon",
    values_to = "abundance",
  ) %>%
  select(id, source, taxon, abundance)

data_pivoted$taxon = factor(data_pivoted$taxon, levels = taxa_order)
data_pivoted_A = data_pivoted[data_pivoted$source == "grey" | data_pivoted$source == "lightblue" | data_pivoted$source == "white" ,]
data_pivoted_B = data_pivoted[data_pivoted$source == "B",]

plot_A = ggplot(data_pivoted_A, aes(x = id, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample", y = "Taxonomy of Assigned Reads", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                                             rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                                                                                                   "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                                                                                                   "WOR-3", "Other Taxa")), values = 
                                                                                             rev(c( color_euk,  color_alpha,  color_firm,  color_bacteroid,
                                                                                                    color_desulfo,   
                                                                                                    color_fibro, color_fuso, color_gamma, color_spiro,
                                                                                                    color_pvc, color_weird, color_other))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, color = c("red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "white", "blue", "blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")))  +  guides(fill = guide_legend(reverse = TRUE))



plot_B = ggplot(data_pivoted_B, aes(x = id, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample", y = "Taxonomic Assessment of Reads", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                                            rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                                                                                                   "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                                                                                                   "WOR-3", "Other Taxa")), values = 
                                                                                             rev(c( color_euk,  color_alpha,  color_firm,  color_bacteroid,
                                                                                                    color_desulfo,   
                                                                                                    color_fibro, color_fuso, color_gamma, color_spiro,
                                                                                                    color_pvc, color_weird, color_other))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +  guides(fill = guide_legend(reverse = TRUE))

svg("figures/read_taxonomy.svg", width = 7, height = 5)
plot_A
dev.off()

svg("read_taxonomy.svg", width = 7, height = 5)
# Create the barplot
theme_set(theme_minimal())
plot_grid(plot_grid(
  plot_A + theme(legend.position = "none"), plot_B + theme(legend.position = "none"),
  ncol = 1, labels = c("A", "B"), align = "hv", axis = "b"), get_legend(plot_A), rel_widths = c(0.25, 0.1))
dev.off()


########

# Read data from file
data <- subset(read.csv("data/read_taxonomy.tab", stringsAsFactors = FALSE, sep="\t"), select = -c(Other.Taxa))
data$id = factor(data$id, levels = sample_order)

taxa_order = rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                   "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                   "WOR.3"))

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
  color_euk = "#696969"
  color_other = "#D3D3D3"
}

# Pivot data using tidyverse
data_pivoted <- data %>%
  pivot_longer(
    cols = taxa_order,
    names_to = "taxon",
    values_to = "abundance",
  ) %>%
  select(id, source, taxon, abundance)

data_pivoted$taxon = factor(data_pivoted$taxon, levels = taxa_order)
data_pivoted_A = data_pivoted[data_pivoted$source == "A",]
data_pivoted_B = data_pivoted[data_pivoted$source == "B",]


plot_A = ggplot(data_pivoted_A, aes(x = id, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Sample", y = "Rel. Abundance", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                                         rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                                                                                               "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                                                                                               "WOR-3")), values = 
                                                                                         rev(c( color_euk,  color_alpha,  color_firm,  color_bacteroid,
                                                                                                color_desulfo,   
                                                                                                color_fibro, color_fuso, color_gamma, color_spiro,
                                                                                                color_pvc, color_weird))) +
  theme_bw() + scale_y_continuous(labels = scales::percent_format()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +  guides(fill = guide_legend(reverse = TRUE))
plot_A

plot_B = ggplot(data_pivoted_B, aes(x = id, y = abundance, fill = taxon)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Sample", y = "Rel. Abundance", fill = "Read Taxonomy") + scale_fill_manual(labels = 
                                                                                         rev(c("Eukaryota", "Alphaproteobacteria", "Bacillota", "Bacteroidota", "Desulfovibrionales", "Fibrobacterota",
                                                                                               "Fusobacteriota",  "Gammaproteobacteria", "Spirochaetota", "Verrucomicrobiota", 
                                                                                               "WOR-3")), values = 
                                                                                         rev(c( color_euk,  color_alpha,  color_firm,  color_bacteroid,
                                                                                                color_desulfo,   
                                                                                                color_fibro, color_fuso, color_gamma, color_spiro,
                                                                                                color_pvc, color_weird))) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +  guides(fill = guide_legend(reverse = TRUE))

svg("read_taxonomy_noother.svg", width = 7, height = 5)
# Create the barplot
theme_set(theme_minimal())
plot_grid(plot_grid(
  plot_A + theme(legend.position = "none"), plot_B + theme(legend.position = "none"),
  ncol = 1, labels = c("A", "B"), align = "hv", axis = "b"), get_legend(plot_A), rel_widths = c(0.25, 0.1))
dev.off()


####

# Read data from file
data <- read.csv("data/bin_quality.txt", stringsAsFactors = FALSE, sep="\t")
sample_order = c("R1_6", "R1_12", "R1_24", "R1_25", "R1_28", "R1_29", "R1_38", "R1_43", "R2_8", "R2_24", "R2_26",
                 "F5GI_2", "F5GI_3", "F5HG_2", "F5HG_3", "F6GI_3", "F6GI_4", "F6HG_2", "F6HG_3",
                 "F7GI_2", "F7GI_3", "F7HG_2", "F7HG_3", "F8GI_2", "F8GI_3", "F8HG_2", "F8HG_3")
data$id = factor(data$id, levels = sample_order)

data_long = data %>% 
  pivot_longer(
    cols = c("Medium.Quality", "High.Quality"),
    names_to = "quality",
    values_to = "bins"
  )

svg("bin_count.svg", width = 8, height = 4.5)
ggplot(data_long, aes(fill=quality, y=bins, x=id)) + 
  geom_bar(position="stack", stat="identity") + labs(x = "Sample", y = "Number of MAGs", fill = "MAG Quality") +
  scale_fill_manual(labels = 
                      c("High Quality", "Medium Quality"), values = c("Orange", "skyBlue")
                      ) +
  
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                     panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



#########

# Read data from file
data_scfa <- subset(read.csv("data/scfa_breakdown.tab", stringsAsFactors = FALSE, sep="\t"))
data_scfa = melt(data_scfa, id = c("Taxa"))
data_scfa$Taxa = factor(data_scfa$Taxa, levels = c("Bacillota", "Bacteroidota", "Gammaproteobacteria", "Desulfovibrionales", "Verrucomicrobiota"))
data_scfa$value = data_scfa$value * 100

#red", "darkturquoise", "darkgreen"
colours = c("red", "red", "red", "red", "red", "red", "red", "darkgreen", "darkgreen","darkgreen", "darkturquoise", "darkturquoise","darkturquoise")

scfa_plot = ggplot(data_scfa, aes(x = variable, y = Taxa)) + 
  geom_point(aes(size = value, fill = variable), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 75), range = c(1, 10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "MAGs with Pathway (%)", fill = "")  +theme_light() + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 11, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 10, face = "bold"), 
        panel.background = element_blank(),  
        legend.position = "right") +  
  scale_fill_manual(values = colours, guide = FALSE) +  
  scale_y_discrete(limits = rev,
                   labels=c(
                     "Bacillota" = "Bacillota (78)",
                     "Bacteroidota" = "Bacteroidota (72)",
                     "Gammaproteobacteria" = "Gammaproteobacteria (31)",
                     "Desulfovibrionales" = "Desulfovibrionales (13)",
                     "Verrucomicrobiota" = "Verrucomicrobiota (6)"
                   )) + 
  
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_x_discrete(labels=c("Pyruvate.to.Acetate.Formate..PFL." = "Pyruvate to Acetate (PFL)",
                            "Pyruvate.to.acetate..PFOR." = "Pyruvate to Acetate (PFOR)",
                            "Ethanolamine.utilization" = "Ethanolamine utilization",
                            "CO2.to.Acetate" = expression("CO" ["2"]*" to Acetate"),
                            "Choline.utilization" = "Choline utilization",
                            "Glycine.to.Acetate" = "Glycine to Acetate",
                            "R.pyruvate.to.R.acetate" = "R-pyruvate to R-acetate (porA)",
                            "Succinate.to.Propionate" = "Succinate to Propanoate",
                            "Acrylate.to.Propionate" = "Acrylate to Propanoate",
                            "Threonine.to.Propionate" = "Threonine to Propanoate",
                            "Acetate.to.Butyrate" = "Acetate to Butyrate",
                            "Glutamate.to.Butyrate" = "Glutamate to Butyrate",
                            "Lysine.Degradation" = "Lysine degradation"))



svg("figures/scfa_table.svg", width = 9, height = 5)
scfa_plot
dev.off()
