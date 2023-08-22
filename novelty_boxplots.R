# Section: load packages
{
  list.of.packages = c('RColorBrewer', 'ggplot2', 'ggbeeswarm', 'cowplot', 'stringr',
                       'dplyr', 'tidyr', 'tidyverse')
  
  # Install packages if missing
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  # Load packages
  sapply(list.of.packages, library, character.only = TRUE)
}


# Create a dataframe of all CAZyme and SulfAtlas annotations
df = as.data.frame(read.csv("data/cazy_sulf_nr.tsv", sep = "\t"))
df$annot = as.factor(df$annot)
df$nr_id = as.numeric(df$nr_id)
df$self_id = as.numeric(df$self_id)
df$sulf_id = as.numeric(df$sulf_id)
df$dbCAN_id = as.numeric(df$dbCAN_id)

# Get only the sulfatase annotations
df_sulf = df[str_detect(df$annot, "S1_"), ]

# tabulate the id variable
tab <- table(df_sulf$annot)
# Get the names of the ids that we care about.
# In this case the ids that occur >= 20 times
idx <- names(tab)[tab >= 20]
idx = c("S1_13", "S1_72", "S1_19", "S1_8", "S1_24", "S1_4", "S1_11", "S1_20", "S1_15", "S1_16", "S1_17", "S1_28", "S1_25", "S1_14", "S1_27")

# Only look at the data that we care about
df_sulf2 = df_sulf[df_sulf$annot %in% idx,]

df_sulf_pivot = df_sulf2 %>%
  pivot_longer(nr_id:sulf_id, names_to = "blast", values_to = "id")
df_sulf_pivot$blast =  factor(df_sulf_pivot$blast, levels = c("self_id", "nr_id", "sulf_id"))
df_sulf_pivot$annot = factor(df_sulf_pivot$annot,
                            idx)

# Get counts of how many genes in each section
sulf_summary = df_sulf_pivot %>% group_by(annot) %>% tally()
sulf_summary$n = sulf_summary$n / 3
sulf_summary$blast = "nr_id"


svg("figures/sulf_boxplot.svg",width = 8, height = 4)
p <- ggplot(df_sulf_pivot, aes(annot, id, fill=blast))
p + geom_boxplot() + theme_bw() + scale_fill_discrete(labels=c('Kyphosid Gut', 'Genbank nr', 'SulfAtlas')) + guides(fill= guide_legend(title = "Database"))  +  xlab("Sulfatase Subfamily") + 
  ylab("Percent Identity") + theme(axis.text.x = element_text(size = 11, color = "black", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 11, color = "black")) + ylim(20, 105) + 
  geom_text(data = sulf_summary,
            aes(annot, Inf, label = n), vjust = 1.5)
dev.off()

##### CAZymes
valid_cazy = c("GH16_16", "GH86", "GH16_11", "GH10", "GH82","GH50", "GH117", "GH150", "GH16_12", "PL40", "PL6", "PL38", "GH29", "PL8", "GH95")
cazy_colors = c("red","red","red","red","red","red","red","red","red","red","brown","brown","brown","brown","brown","brown","brown")
df_cazy = df[df$annot %in% valid_cazy,]
df_cazy_pivot = df_cazy %>%
  pivot_longer(c(nr_id, self_id, dbCAN_id), names_to = "blast", values_to = "id")
df_cazy_pivot$blast =  factor(df_cazy_pivot$blast, levels = c("self_id", "nr_id", "dbCAN_id"))
df_cazy_pivot$annot <- factor(df_cazy_pivot$annot,     # Reorder factor levels
                             valid_cazy)

cazy_summary = df_cazy_pivot %>% group_by(annot) %>% tally()
cazy_summary$n = cazy_summary$n / 3
cazy_summary$blast = "nr_id"


svg("figures/cazy_boxplot.svg", width = 8, height = 4)
p <- ggplot(df_cazy_pivot, aes(annot, id, fill=blast))
p + geom_boxplot() + theme_bw() + scale_fill_discrete(labels=c('Kyphosid Gut', 'Genbank nr', 'CAZy db')) + guides(fill= guide_legend(title = "Database"))  +  xlab("CAZyme Class") + 
  ylab("Percent Identity")+ theme(axis.text.x = element_text(size = 11, color = cazy_colors,angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size = 11, color = "black")) + ylim(20, 105) +
  geom_text(data = cazy_summary,
            aes(annot, Inf, label = n), vjust = 1.5)
dev.off()


