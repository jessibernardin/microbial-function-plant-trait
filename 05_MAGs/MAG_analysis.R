#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : June 25, 2024 ####
#### Metagenome Assembled Genome (MAG) Analysis

# Load and install required packages
packages_to_load <- c(
  "tidyverse", "dplyr", "pheatmap", "qiime2R", "phyloseq", "this.path", "magrittr"
)

for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

setwd(this.path::here())

#### MAPPING TRANSCRIPTS AND FUNCTIONS TO MAGS ####
#read in csv blast results of all transcripts blasted against drep_mag blast database
mag_blast <- read.delim("blast_output_formatted.csv", header = FALSE)
colnames(mag_blast) <- c("contig", "contig_bin", "pident", "bitscore", "evalue")
mag_qual <- read.csv("dRep_checkm2_quality_report.csv", header=TRUE)

mag_tax <- read.csv("dRep_mag_tax.csv", header=TRUE)
new_bin <- mag_tax[,1]
mag_tax <- mag_tax[,-1]
names(mag_tax)[1] <- "Taxon"
names(mag_tax)[2] <- "Feature.ID"

mag_tax_p <- parse_taxonomy(mag_tax)
mag_tax_p$bin <- rownames(mag_tax_p)
mag_tax_p <- cbind(new_bin, mag_tax_p)
mag_tax_filt<- mag_tax_p[,c(1,8,9)]
split_column <- strsplit(mag_blast$contig_bin, "_", fixed = TRUE)

# Create new columns for contig and bin
mag_blast$gene <- sapply(split_column, function(x) paste(x[1:2], collapse = "_"))
mag_blast$bin <- sapply(split_column, function(x) paste(x[-(1:2)], collapse = "_"))

#functional annotation
kegg.f <- read.csv("kegg.f.csv", header=TRUE)
kegg.f <- kegg.f[,-1]
#merge with KEGG_KO
mag_fun <- left_join(mag_blast, kegg.f, by = c("contig" = "gene_name"))
#260,095

mag_fun_ko_def <- mag_fun[,c(7,8,12)]

filtered_df <- mag_fun_ko_def %>% filter(grepl("aminopeptidase|chit|acetylhexosaminidase", KO_definition, ignore.case = TRUE))

filtered_df <- merge(filtered_df, mag_tax_p, by.x="bin", by.y="bin")
filtered_df$mag_species <- paste(filtered_df$new_bin, filtered_df$Species, sep = "-")
filtered_df <- filtered_df %>%
  mutate(broad_enzyme = case_when(
    str_detect(KO_definition, regex("amino", ignore_case = TRUE)) ~ "Protease",
    str_detect(KO_definition, regex("chit|acetylhexosaminidase", ignore_case = TRUE)) ~ "Chitinase",
    TRUE ~ "Other"  # Handle cases where none of the patterns match
  ))

detach("package:plyr", unload=TRUE)
filtered_df_new <- filtered_df %>%
  group_by(mag_species, KO_definition) %>%
  mutate(value = n()) %>%
  ungroup()

colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
             "#F6AE2D","#86BBD8")
colours2 = c( "#F7EE7F","#679289")

filtered_df_new$value <- as.integer(filtered_df_new$value)



filtered_df_new <- filtered_df_new[order(filtered_df_new$value, decreasing = TRUE),]  

filtered_df_new$broad_enzyme <- factor(filtered_df_new$broad_enzyme, levels = c("Chitinase", "Protease"))

filtered_df_new <- filtered_df_new %>%
  arrange(broad_enzyme)

filtered_df_new$KO_definition <- as.factor(filtered_df_new$KO_definition)
ggplot(filtered_df_new, aes(y = reorder(mag_species, value), x = factor(KO_definition, levels=unique(KO_definition)), size = value, fill = KO_definition)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 10), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "none")

filtered_df_2 <- filtered_df_new

filtered_df_2 <- filtered_df_2 %>%
  group_by(mag_species, broad_enzyme) %>%
  mutate(value_be = n()) %>%
  ungroup()

ggplot(filtered_df_2, aes(x = reorder(mag_species, -value_be), y = broad_enzyme, size = value_be, fill = broad_enzyme)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 17), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, angle = 45, vjust = 1.01, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "right") +  
  scale_fill_manual(values = colours2, guide = FALSE)+coord_flip()


#### MAG PHYLOGENETICS ####
library(ggtree)
library(ape)

# Read the phylogenetic tree
tree <- read.tree("aligned_marker_genes.fasta.treefile")
tree$tip.label <- gsub("_modified_drep", "", tree$tip.label)

# Plot the tree
ggtree(tree) + geom_tiplab()


# Extract tip labels (order of MAGs)
mag_order <- tree$tip.label

# Assuming counts_long contains the relative abundance data with a column for MAG IDs
# Arrange the data according to the phylogenetic order
filtered_df_new$bin <- factor(filtered_df_new$bin, levels = mag_order)

#Fig8
ggplot(filtered_df_new, aes(y = mag_species, x = factor(KO_definition, levels=unique(KO_definition)), size = value, fill = KO_definition)) + 
  geom_point(alpha = 0.75, shape = 21, stroke = NA) + 
  scale_size_continuous(limits = c(1, 300), range = c(1, 10), breaks = c(10, 50, 100, 300)) + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "none")


#### MAG relative abundance ####

# Load the combined counts
counts <- read.table("mapped_read_counts.txt", header=TRUE, sep="\t")
counts<- counts[-1,]
counts <- as.data.frame(counts)
counts$MappedReadCount <- as.integer(counts$MappedReadCount)

split_sample <- function(sample_column, delimiter, n) {
  parts <- str_split_fixed(sample_column, delimiter, n)
  prefix <- apply(parts[, 1:(n-1), drop = FALSE], 1, paste, collapse = delimiter)
  suffix <- parts[, n]
  data.frame(Sample_ID = prefix, bin = suffix)
}

df_split <- split_sample(counts$Sample, "_", 5)
df_combined <- cbind(counts, df_split)
df_combined <- df_combined[,-1]

df_aggregated <- df_combined %>%
  group_by(Sample_ID, bin) %>%
  summarise(TotalMappedReads = sum(MappedReadCount)) %>%
  ungroup() %>%
  mutate(RelativeAbundance = TotalMappedReads / sum(TotalMappedReads))


separated_df <- df_combined %>%
  separate(Sample_ID, into = c("sample", "plant", "treatment", "week"), sep = "_")

separated_df$treatment <- ifelse(grepl("M01", separated_df$treatment) , "CommA",
                                 ifelse(grepl("M06", separated_df$treatment) , "CommB",
                                        ifelse(grepl("M09", separated_df$treatment) , "CommC",NA)))


separated_df <- separated_df %>%
  filter(!(bin %in% c("concoct.74", "maxbin.41")))


df_aggregated2 <- separated_df %>%
  group_by(treatment, bin) %>%
  summarise(TotalMappedReads = sum(MappedReadCount)) %>%
  ungroup() %>%
  group_by(treatment) %>%
  mutate(RelativeAbundance = TotalMappedReads / sum(TotalMappedReads))


df_aggregated3 <- left_join(df_aggregated2, mag_tax_filt, by="bin")
df_aggregated3 <- df_aggregated3 %>%
  filter(!is.na(new_bin))



ggplot(df_aggregated3, aes(x =treatment , y =Species , size = RelativeAbundance, fill = Species)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous() + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "none")



df_aggregated3$bin <- factor(df_aggregated3$bin, levels = mag_order)


#FIG8
ggplot(df_aggregated3, aes(x = treatment, y = new_bin, fill = RelativeAbundance)) +
  geom_tile(color = "white") +  # Use geom_tile() for heatmap
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0)+
  labs(x = "", y = "", fill = "Relative Abundance") +  # Label adjustments
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    legend.position = "right"  # Position legend
  )



#### redo relative abundance based on metagenomic mapping instread of metatranscr.


# Load the combined counts
counts <- read.table("megaG_mag_abund/mapped_read_counts.txt", header=TRUE, sep="\t")
counts<- counts[-1,]
counts <- as.data.frame(counts)
counts$MappedReadCount <- as.integer(counts$MappedReadCount)

split_sample <- function(sample_column, delimiter, n) {
  parts <- str_split_fixed(sample_column, delimiter, n)
  prefix <- apply(parts[, 1:(n-1), drop = FALSE], 1, paste, collapse = delimiter)
  suffix <- parts[, n]
  data.frame(Sample_ID = prefix, bin = suffix)
}

df_split <- split_sample(counts$Sample, "_", 5)
df_combined <- cbind(counts, df_split)
df_combined <- df_combined[,-1]

df_aggregated <- df_combined %>%
  group_by(Sample_ID, bin) %>%
  summarise(TotalMappedReads = sum(MappedReadCount)) %>%
  ungroup() %>%
  mutate(RelativeAbundance = TotalMappedReads / sum(TotalMappedReads))

 
separated_df <- df_combined %>%
  separate(Sample_ID, into = c("sample", "plant", "treatment", "week"), sep = "_")

separated_df$treatment <- ifelse(grepl("M01", separated_df$treatment) , "CommA",
                                 ifelse(grepl("M06", separated_df$treatment) , "CommB",
                                        ifelse(grepl("M09", separated_df$treatment) , "CommC",NA)))

separated_df$bin <- gsub("_modified_drep_counts", "", separated_df$bin)
separated_df$bin <- gsub("_modified_drep.bam.tmp.0000_counts", "", separated_df$bin)
separated_df$bin <- gsub("_modified_drep.bam.tmp.0001_counts", "", separated_df$bin)
separated_df <- separated_df[-614,]


separated_df <- separated_df %>%
  filter(!(bin %in% c("concoct.74", "maxbin.41")))


df_aggregated2 <- separated_df %>%
  group_by(treatment, bin) %>%
  summarise(TotalMappedReads = sum(MappedReadCount)) %>%
  ungroup() %>%
  group_by(treatment) %>%
  mutate(RelativeAbundance = TotalMappedReads / sum(TotalMappedReads))

ggplot(df_aggregated2, aes(x = treatment, y = RelativeAbundance, fill = bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Treatment", y = "Relative Abundance", title = "Relative Abundance of MAGs by Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
colours3 = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
              "#F6AE2D","#86BBD8","#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
              "#F6AE2D","#86BBD8","#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
              "#F6AE2D","#86BBD8","#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
              "#F6AE2D","#86BBD8")


ggplot(df_aggregated2, aes(x = treatment, y = RelativeAbundance, fill = bin)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Treatment", y = "Relative Abundance", title = "Relative Abundance of MAGs by Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=colours3)


df_aggregated3 <- left_join(df_aggregated2, mag_tax_filt, by="bin")
df_aggregated3 <- df_aggregated3 %>%
  filter(!is.na(new_bin))
ggplot(df_aggregated3, aes(x = treatment, y = RelativeAbundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Treatment", y = "Relative Abundance", title = "Relative Abundance of MAGs by Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_fill_manual(values=colours3)



ggplot(df_aggregated3, aes(x =treatment , y =Species , size = RelativeAbundance, fill = Species)) + 
  geom_point(alpha = 0.75, shape = 21) + 
  scale_size_continuous() + 
  labs(x = "", y = "", size = "# Unique Contigs", fill = "") +theme_classic()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, vjust = 1.01, hjust = 1, angle=90), 
        axis.text.y = element_text(colour = "black", size = 5), 
        legend.text = element_text(size = 5, colour ="black"), 
        legend.position = "none")



df_aggregated3$bin <- factor(df_aggregated3$bin, levels = mag_order)



ggplot(df_aggregated3, aes(x = treatment, y = new_bin, fill = RelativeAbundance)) +
  geom_tile(color = "white") +  # Use geom_tile() for heatmap
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0)+
  labs(x = "", y = "", fill = "Relative Abundance") +  # Label adjustments
  theme_minimal() +  # Minimal theme
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # Rotate x-axis labels
    legend.position = "right"  # Position legend
  )


















