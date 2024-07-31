#### Microbial function mediates leaf traits in a pitcher plant model system ####
#### Authors: Jessica R. Bernardin, Erica B. Young, Leonora S. Bittleston ####
#### last update : June 21, 2024 ####
#### Metatranscriptomic Functional Analysis

#### LOAD PACKAGES AND SET WD ####
packages_to_load <- c(
  "tidyr", "dplyr", "compositions", "ggpubr", "variancePartition", "edgeR", "BiocParallel",
  "stringr", "pheatmap", "qiime2R","this.path", "limma", "ggrepel", "vegan", "VennDiagram", "ggforce",
  "tidyverse", "knitr"
)

for (i in packages_to_load) { #Installs packages if not yet installed
  if (!require(i, character.only = TRUE)) install.packages(i)
}

#BiocManager::install("variancePartition")
#BiocManager::install("edgeR")

setwd(this.path::here())

#### CLUSTER REDUNDANT TRANSCRIPTS USING CDHIT MAP ####
#https://rpubs.com/rmurdoch/cdhit_to_mapping_file
clstr <- read.csv("metaTfilt2_faa_cdhit.faa.clstr", sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
clstr2 <- clstr
n = nrow(clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
  if (numbers_only(clstr2[row,1]) == TRUE) {
    clstr2[row,1] <- x}
  else {NULL}
  x <- clstr2[row,1]
}

clstr4 <- clstr2[-which(clstr2$V2 == ""), ]

clstr5 <- clstr4
clstr5[] <- lapply(clstr5, gsub, pattern='>', replacement='')
clstr5.2 <- data.frame(str_split_fixed(clstr5$V2, "aa, ", 2))
clstr5.3 <- data.frame(str_split_fixed(clstr5.2$X2, "... ", 2))
clstr6 <- cbind(clstr5[1],clstr5.2[1],clstr5.3[1:2])
colnames(clstr6) <- c("cluster","aa","gene","stat")
clstr6 <- clstr6 %>%
  mutate(rep = ifelse(stat == "*", TRUE, FALSE))

write_csv(clstr6, "cdhit_map.csv")

#### MAKE COUNTS AND TPM TABLES FROM KALLISTO OUTPUTS ####

##COUNT
files <- list.files(path = "kallisto/", pattern = "\\.tsv$", full.names = TRUE)
df_list <- list()
for (file in files) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_list[[file]] <- data[[4]]
}
counts <- do.call(cbind, df_list)
data_first_file <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_id <- data_first_file[[1]]
length <- data_first_file[[2]]
counts <- cbind(target_id, length, counts)
counts <- as.data.frame(counts)
rownames(counts) <- counts$target_id
counts[,-1] <- apply(counts[,-1], 2, function(x) as.numeric(as.character(x)))
counts <- counts[,-c(1:2)]
new_column_names <- gsub("^kallisto//|\\-aligned\\.abundance\\.tsv$", "", colnames(counts))
colnames(counts) <- new_column_names
sum(rowSums(counts) == 0)#8898
write.csv(counts, "counts.csv")
counts <- read.csv("counts.csv", header=TRUE)


##TPM
for (file in files) {
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df_list[[file]] <- data[[5]]
}
tpm <- do.call(cbind, df_list)
data_first_file <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
target_id <- data_first_file[[1]]
length <- data_first_file[[2]]
tpm <- cbind(target_id, length, tpm)
tpm <- as.data.frame(tpm)
rownames(tpm) <- tpm$target_id
tpm[,-1] <- apply(tpm[,-1], 2, function(x) as.numeric(as.character(x)))
tpm <- tpm[,-c(1:2)]
new_column_names <- gsub("^kallisto//|\\-aligned\\.abundance\\.tsv$", "", colnames(tpm))
colnames(tpm) <- new_column_names
sum(rowSums(tpm) == 0)#8898
write.csv(tpm, "tpm.csv")
tpm <- read.csv("tpm.csv", header=TRUE)
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]

###collapse by cluster
cluster <- read.csv("cdhit_map.csv", header=TRUE)
tpm$gene <- rownames(tpm)
tpm_cluster <- left_join(cluster, tpm, by="gene")
names(tpm_cluster)
rna_columns <- grep("RNA", names(tpm_cluster), value = TRUE)
tpm_cluster_sum <- tpm_cluster %>%
  group_by(cluster) %>%
  summarize(
    gene = gene[rep == TRUE][1],
    across(starts_with(rna_columns), sum, na.rm = TRUE),
    .groups = 'drop')

colSums(tpm_cluster_sum[,c(3:29)])
dim(tpm_cluster_sum)#225735      29
sum(rowSums(tpm_cluster_sum[,3:29]) == 0)#7398

tpm_cluster_sum <- as.data.frame(tpm_cluster_sum)
rownames(tpm_cluster_sum) <- tpm_cluster_sum$gene

tpm_cluster_sum <- tpm_cluster_sum[,-c(1,2)]
tpm_cluster_sum <- round(tpm_cluster_sum, digits = 0)
str(tpm_cluster)
write.csv(tpm_cluster_sum, "tpm_cluster_sum.csv")
tpm_cluster_sum <- read.csv("tpm_cluster_sum.csv", header=TRUE)
rownames(tpm_cluster_sum) <- tpm_cluster_sum[,1]
tpm_cluster_sum <- tpm_cluster_sum[,-1]

#COUNTS
counts$gene <- rownames(counts)
count_cluster <- left_join(cluster, counts, by="gene")
names(count_cluster)
rna_columns <- grep("RNA", names(count_cluster), value = TRUE)
counts_cluster_sum <- count_cluster %>%
  group_by(cluster) %>%
  summarize(
    gene = gene[rep == TRUE][1],
    across(starts_with(rna_columns), sum, na.rm = TRUE),
    .groups = 'drop')

colSums(counts_cluster_sum[,c(3:29)])
dim(counts_cluster_sum)#225,735     29
sum(rowSums(counts_cluster_sum[,3:29]) == 0)#7398
counts_cluster_sum <- as.data.frame(counts_cluster_sum)
rownames(counts_cluster_sum) <- counts_cluster_sum$gene
counts_cluster_sum <- counts_cluster_sum[,-c(1,2)]
counts_cluster_sum <- round(counts_cluster_sum, digits = 0)
write.csv(counts_cluster_sum, "counts_cluster.csv")
counts_cluster_sum <- read.csv("counts_cluster.csv", header=TRUE)
rownames(counts_cluster_sum) <- counts_cluster_sum[,1]
counts_cluster_sum <- counts_cluster_sum[,-1]



#read in metadata
meta_all <- read.csv("Exp_1_function_metadata_2021.csv", header=TRUE)
meta_rna <- meta_all %>% filter(nucleic_acid == "RNA")
meta_rna <- meta_rna[,-1]
rownames(meta_rna) <- meta_rna[,1]
meta_rna$treatment <- as.factor(meta_rna$treatment)

#removing samples with low reads from the DREAM analysis
meta_rna2 <- meta_rna %>%
  filter(!sample_id %in% c("RNA5_P3_M01_WK4", "RNA9_P13_M01_WK9", "RNA10_P5_M06_WK1", "RNA18_P45_M06_WK9"))
rownames(meta_rna2) <- meta_rna2[,1]

#### READ IN AND FORMAT FUNCTIONAL ANNOTATION FROM KOFAMSCAN ####
ko_tsv <- read.table("metaTfilt2.kofamscan_detailed.tsv", sep = "\t", header = FALSE) #5,344,716 genes/kos
kegg.f <- ko_tsv
kegg.f<-subset(kegg.f, V6 < 0.001)[,c(2:7)] #e value less than .001, 3,342,730

#find the KO with the best match for each transcript (using bit score)
kegg.f <- kegg.f %>%
  group_by(V2) %>%
  filter(V5 == max(V5)) %>%
  ungroup()
colnames(kegg.f) <- c("gene_name", "KO", "thrshld", "score", "E_value", "KO_definition")#174,370 unique transcripts mapped to one function


write.csv(kegg.f, "../05_MAGs/kegg.f.csv")
#read in KEGG orthology from website, has pathway info for the KEGGs
KEGG_hier <- read.csv("KO_Orthology_ko00001.csv", header=FALSE) #https://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
names(KEGG_hier) <- c("L1", "L2", "L3", "L4")
KEGG_hier$KO <- word(KEGG_hier$L4, 1) #61,302 KOs
KEGG_hier$KO <- as.factor(KEGG_hier$KO)
kegg.f$KO <- as.factor(kegg.f$KO)
KEGG_KO <- left_join(kegg.f, KEGG_hier, by="KO") #352378

#reorder the rows and columns
matched_cols <- colnames(tpm_cluster_sum)[match(rownames(meta_rna2), colnames(tpm_cluster_sum))]
tpm_cluster_sum <- tpm_cluster_sum[, matched_cols]
colnames(tpm_cluster_sum) == rownames(meta_rna2)

matched_cols <- colnames(counts_cluster_sum)[match(rownames(meta_rna2), colnames(counts_cluster_sum))]
counts_cluster_sum <- counts_cluster_sum[, matched_cols]
colnames(counts_cluster_sum) == rownames(meta_rna2)
#remove genes whose rowsums are 0 for counts
sum(rowSums(counts_cluster_sum) == 0)
#7418
dim(counts_cluster_sum)
#225735
counts_cluster_sum_filt<- counts_cluster_sum[rowSums(counts_cluster_sum) != 0, ]
dim(counts_cluster_sum_filt)
#218317

#remove genes whose rowsums are 0 for tpm
sum(rowSums(tpm_cluster_sum) == 0)
#7418
dim(tpm_cluster_sum)
#225735
tpm_cluster_sum_filt<- tpm_cluster_sum[rowSums(tpm_cluster_sum) != 0, ]
dim(tpm_cluster_sum_filt)
#218317

#### FILTER GENES OUT THAT HAVE NO FUNCTIONAL ANNOTATION

##COUNTS
counts_fun <- counts_cluster_sum_filt
counts_fun$gene_name <- rownames(counts_fun)
kegg.f <- kegg.f %>%
  distinct(gene_name, .keep_all = TRUE)
counts_fun_comb <- dplyr::left_join(counts_fun, kegg.f, by="gene_name")
# 218317 genes and rows combined with functions

#remove rows with no function
counts_fun_filt <- counts_fun_comb %>% filter(!is.na(KO))
# 162624 remain after filtering, meaning that xxxx genes were missing functional annotation

#subset the df to just those columns we need
counts_fun_meta <- counts_fun_filt[,c(24:29)]
counts_fun_subset <- counts_fun_filt[,c(1:24)]
rownames(counts_fun_subset) <- counts_fun_subset$gene_name
counts_fun_subset <- counts_fun_subset[,-24]

##TPM
tpm_fun <- tpm_cluster_sum_filt
tpm_fun$gene_name <- rownames(tpm_fun)
tpm_fun_comb <- left_join(tpm_fun, kegg.f, by="gene_name")
# 218317 genes and rows combined with functions

#remove rows with no function
tpm_fun_filt <- tpm_fun_comb %>% filter(!is.na(KO))
# 162623 remain after filtering, meaning that xxxx genes were missing functional annotation

#subset the df to just those columns we need
tpm_fun_meta <- tpm_fun_filt[,c(24:29)]
tpm_fun_subset <- tpm_fun_filt[,c(1:24)]
rownames(tpm_fun_subset) <- tpm_fun_subset$gene_name
tpm_fun_subset <- tpm_fun_subset[,-24]


#match up with the metadata
matched_cols <- colnames(counts_fun_subset)[match(rownames(meta_rna2), colnames(counts_fun_subset))]
counts_fun_subset <- counts_fun_subset[, matched_cols]
colnames(counts_fun_subset) == rownames(meta_rna2)

#match up with the metadata
matched_cols <- colnames(tpm_fun_subset)[match(rownames(meta_rna2), colnames(tpm_fun_subset))]
tpm_fun_subset <- tpm_fun_subset[, matched_cols]
colnames(tpm_fun_subset) == rownames(meta_rna2)

#### RUN DIFFERENTIAL EXPRESSION MODEL *DREAM* ####
dge2 <- DGEList(counts_cluster_sum_filt) ##the matrix with gene names as row names and counts in columns, 218317 (was 83977) genes
norm.mat2 <- calcNormFactors(dge2)

dge <- DGEList(counts_fun_subset) ##the matrix with gene names as row names and counts in columns, 162,624, FILTERED OUT GENES WITH NO FUNCTIONAL ANNOTATION
norm.mat <- calcNormFactors(dge)

# Specify parallel processing parameters
param <- SnowParam(4, "SOCK", progressbar=TRUE)
param <- SnowParam(workers = parallel::detectCores() - 1, type = "SOCK", progressbar = TRUE)

form <- ~ 0 + treatment + (1|plant_number)

# estimate weights using linear mixed model of dream (no intercept baseline)
vobjDream2 <- voomWithDreamWeights(norm.mat2, form, meta_rna2, BPPARAM=param )

vobjDream <- voomWithDreamWeights(norm.mat, form, meta_rna2, BPPARAM=param )

#make contrasts
L1 <- makeContrastsDream(form, meta_rna2, 
                         contrasts = c(TestCommA = "treatmentCommA - (treatmentCommB + treatmentCommC)/2",
                                       TestCommB = "treatmentCommB - (treatmentCommA + treatmentCommC)/2",
                                       TestCommC = "treatmentCommC - (treatmentCommB + treatmentCommA)/2"))

#model 218317 genes (these were clustered (summed) based on redundancy)
fitm1 <- dream(vobjDream2, form, meta_rna2, L=L1, BPPARAM=param)
fitm1 <- eBayes(fitm1)
saveRDS(fitm1, "dream_model_cluster_06_19_24.RDS")
fitm1 <- readRDS("dream_model_cluster_06_19_24.RDS")
head(fitm1$coef, 3)


#model 162,624 genes (these were clustered based on redundancy and genes with no KO removed
fitm2 <- dream(vobjDream, form, meta_rna2, L=L1, BPPARAM=param)
attr(fitm2, 'errors')
#no convergence for k141_17275_2, k141_118850_2, k141_93065_2

fitm2 <- eBayes(fitm2)
saveRDS(fitm2, "dream_model_NA_FUN_REMOVED_06_21_24.RDS")
fitm2 <- readRDS("dream_model_NA_FUN_REMOVED_06_21_24.RDS")
head(fitm2$coef, 3)

# comparing the effect to the mean of the others
topA2 <- topTable(fitm2, coef = "TestCommA", n=Inf)
topB2 <- topTable(fitm2, coef = "TestCommB", n=Inf)
topC2 <- topTable(fitm2, coef = "TestCommC", n=Inf)

topA <- topA2
topB <- topB2
topC <- topC2

topA$treatment_comp <- "CommA"
topB$treatment_comp <- "CommB"
topC$treatment_comp <- "CommC"

topA$gene_name <- rownames(topA)
topB$gene_name <- rownames(topB)
topC$gene_name <- rownames(topC)

top_all <- rbind(topA, topB, topC)
top_all <- as.data.frame(top_all)

tab1_filtered2 <- topA2[topA2$adj.P.Val < 0.05, ] 
tab2_filtered2 <- topB2[topB2$adj.P.Val < 0.05, ] 
tab3_filtered2 <- topC2[topC2$adj.P.Val < 0.05, ] 

tab1_filtered2$treatment_comp <- "CommA"
tab2_filtered2$treatment_comp <- "CommB"
tab3_filtered2$treatment_comp <- "CommC"

tab1_filtered2$gene_name <- rownames(tab1_filtered2)
tab2_filtered2$gene_name <- rownames(tab2_filtered2)
tab3_filtered2$gene_name <- rownames(tab3_filtered2)

countsCommA <- table(tab1_filtered2$logFC > 0)
#FALSE  TRUE 
#10885  5767 



countsCommB <- table(tab2_filtered2$logFC > 0)
#FALSE  TRUE 
#4278  2124 



countsCommC <- table(tab3_filtered2$logFC > 0)
#FALSE  TRUE 
#6257  5851


tab_filtered2 <- rbind(tab1_filtered2, tab2_filtered2, tab3_filtered2)#110,315 genes, 465 functions, 35,162




#### COMBINING ABUNDANCE WITH LFC DATAFRAMES ####
#match sig genes with TPM and ko
tpm_fun_subset$gene_name <- rownames(tpm_fun_subset)
counts_fun_subset$gene_name <- rownames(counts_fun_subset)


#### Functional Differences ####
unique_genes <- unique(tab_filtered2$gene_name) #24357 unique DE genes

# Filter transcript abundance table based on unique genes
tpm_gene_sig <- tpm_fun_subset[tpm_fun_subset$gene_name %in% unique_genes, ]

#merge with KO functions
tpm_gene_sig_KO <- left_join(tpm_gene_sig, KEGG_KO , by="gene_name")

melted_sig_genes <- tpm_gene_sig_KO %>%
  gather(key = Sample, value = TPM, starts_with("RNA"))

separated_df <- melted_sig_genes %>%
  separate(Sample, into = c("sample", "plant", "treatment", "week"), sep = "_")

separated_df$treatment <- ifelse(grepl("M01", separated_df$treatment) , "CommA",
                                 ifelse(grepl("M06", separated_df$treatment) , "CommB",
                                        ifelse(grepl("M09", separated_df$treatment) , "CommC",NA)))

separated_df$week <- ifelse(grepl("1", separated_df$week) , "0",
                            ifelse(grepl("4", separated_df$week) , "3",
                                   ifelse(grepl("9", separated_df$week) , "8",NA)))

separated_df <- separated_df %>% mutate(Day =
                                          case_when(week <= 0 ~ 1, 
                                                    week <= 3 ~ 22,
                                                    week >= 8 ~ 55))

separated_df <- as.data.frame(separated_df)

sum_tpm_by_treatment_day <- separated_df %>%
  group_by( Day, plant) %>%
  summarize(mean_TPM = mean(TPM), .groups = 'drop')

# Use this function to bin KOs into functional categories
find_category <- function(description) {
  if (grepl("dipeptidase", description, ignore.case = TRUE)) {
    return("Dipeptidase")
  } else if (grepl("aminopeptidase", description, ignore.case = TRUE)) {
    return("Aminopeptidase")
  } else if (grepl("endopeptidase", description, ignore.case = TRUE)) {
    return("Endopeptidase")
  } else if (grepl("peptidase", description, ignore.case = TRUE)) {
    return("Peptidase")
  } else if (grepl("protease", description, ignore.case = TRUE)) {
    return("Protease")
  } else if (grepl("chitinase", description, ignore.case = TRUE)) {
    return("Endochitinase")
  } else if (grepl("endo-chitode", description, ignore.case = TRUE)) {
    return("Endo-chitodextinase")
  } else if (grepl("acetylhexosaminidase", description, ignore.case = TRUE)) {
    return("Exochitinase")
  } else if (grepl("glycosidase|glucosidase", description, ignore.case = TRUE)) {
    return("Glucosidases")
  } else {
    return("Other")  # If none of the terms are found
  }
}

# Apply the function to create the new column
separated_df2 <- mutate(separated_df, KO_category = sapply(KO_definition, find_category))
KEGG_KO2 <- mutate(KEGG_KO, KO_category = sapply(KO_definition, find_category))

#### CHITINASE ####
filtered_chit <- KEGG_KO2 %>%
  filter(str_detect(KO_definition, regex("chit|acetylhexosaminidase", ignore_case = TRUE)))

##chitinase tpm
chit2_filtall <- subset(separated_df2, separated_df2$gene_name %in% filtered_chit$gene_name)
chit2_filtall$newname1 <- paste(chit2_filtall$gene_name, " - ", chit2_filtall$KO_definition)
chit2_filtall <- chit2_filtall[,c(1,2,6,13,15,17,18)]

chit2_filtall_combo <- chit2_filtall %>% 
  group_by(KO_category, newname1, treatment) %>% 
  summarize(tpm_sum=sum(TPM)) %>% 
  ungroup()


unique_genes <- unique(chit2_filtall_combo$newname1)
commB_summary <- aggregate(tpm_sum ~ newname1, data = chit2_filtall_combo[chit2_filtall_combo$treatment == "CommB", ], FUN = sum)  
sorted_genes <- commB_summary[order(commB_summary$tpm_sum, decreasing = FALSE), "newname1"]
genes_not_in_CommB <- unique_genes[!(unique_genes %in% sorted_genes)]
sorted_genes_all <- c(sorted_genes, genes_not_in_CommB)

chit2_filtall_combo$newname1 <- factor(chit2_filtall_combo$newname1, levels = sorted_genes_all)

c2 <- ggplot(chit2_filtall_combo,aes(x = treatment, y = newname1, fill = log2(tpm_sum+0.5))) +
  geom_tile() + 
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),  # Remove gridlines
    panel.background = element_rect(fill = "white", color = NA))+
  facet_col(~KO_category, scales = "free_y", space = "free")+theme(axis.text.y = element_blank())+theme(legend.position = "none")
  
#
#filter the genes to just those associated with enzyme kos then plot lfc as heatmap
chit2_filt1 <- subset(tab1_filtered2, tab1_filtered2$gene_name %in% filtered_chit$gene_name)
chit2_filt2 <- subset(tab2_filtered2, tab2_filtered2$gene_name %in% filtered_chit$gene_name)
chit2_filt3 <- subset(tab3_filtered2, tab3_filtered2$gene_name %in% filtered_chit$gene_name)

chit2_filtall <- rbind(chit2_filt1, chit2_filt2, chit2_filt3)
chit2_filtall_def <- left_join(chit2_filtall, filtered_chit, by="gene_name")
chit2_filtall_def$newname1 <- paste(chit2_filtall_def$gene_name, " - ", chit2_filtall_def$KO_definition)
chit2_filtall_def$newname2 <- paste(chit2_filtall_def$gene_name, " - ", chit2_filtall_def$KO_category)

complete_data <- expand.grid(newname1 = unique(chit2_filtall_def$newname1),
                             treatment_comp = unique(chit2_filtall_def$treatment_comp))

# Merge complete dataset with original data to fill in missing logFC values
chit2_complete <- merge(complete_data, chit2_filtall_def, by = c("newname1", "treatment_comp"), all.x = TRUE)

# Replace missing logFC values with 0
chit2_complete$logFC[is.na(chit2_complete$logFC)] <- 0
chit2_complete$treatment_comp <- factor(chit2_complete$treatment_comp, levels=c("CommA", "CommB", "CommC"))
chit2_complete$newname1 <- factor(chit2_complete$newname1, levels = sorted_genes_all)

c1 <- chit2_complete %>%
  filter(!is.na(score)) %>%ggplot( aes(x = logFC, y = newname1 , color = treatment_comp)) +
  geom_point(position = position_dodge(width = .5), size=4, alpha=.7) +  
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values=c("#014b7a","orange", "#44b7c2")) +
  theme_minimal() +
  labs(x = "Log Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlim(-5, 5)+
  facet_col(~KO_category, scales = "free_y", space = "free")+theme(axis.text.y = element_blank())+theme(legend.position = "none")


#### PROTEASE ####
filtered_prot <- KEGG_KO2 %>%
  filter(str_detect(KO_definition, regex("aminopeptidase|peptid", ignore_case = TRUE)))

# Define a function to search for terms and return the category
find_category2 <- function(description) {
  if (grepl("leucyl aminopeptidase|bacterial leucyl aminopeptidase", description, ignore.case = TRUE)) {
    return("Leucyl Aminopeptidase")
  } else if (grepl("tripeptide aminopeptidase", description, ignore.case = TRUE)) {
    return("Tripeptide Aminopeptidase")
  } else if (grepl("tetrahedral aminopeptidase", description, ignore.case = TRUE)) {
    return("Tetrahedral Aminopeptidase")
  } else if (grepl("methionyl aminopeptidase", description, ignore.case = TRUE)) {
    return("Methionyl Aminopeptidase")
  } else if (grepl("aspartyl aminopeptidase|Xaa-Pro aminopeptidase 2|aminopeptidase YwaD|Xaa-Pro aminopeptidase|PepB aminopeptidase|glutamyl aminopeptidase", description, ignore.case = TRUE)) {
    return("Other")
  } else if (grepl("aminopeptidase N|aminopeptidase|D-aminopeptidase", description, ignore.case = TRUE)) {
    return("Aminopeptidase")
  } else {
    return("Other")  # If none of the terms are found
  }
}

#filter the genes to just those associated with enzyme kos then plot tpm as heatmap
protdf<- separated_df2
protdf <- subset(protdf, KO_category == "Aminopeptidase")
protdf_complete <- mutate(protdf, KO_category2 = sapply(KO_definition, find_category2))

prot2_filtall <- subset(protdf_complete, protdf_complete$gene_name %in% filtered_prot$gene_name)
prot2_filtall$newname1 <- paste(prot2_filtall$gene_name, " - ", prot2_filtall$KO_definition)
prot2_filtall <- prot2_filtall[,c(1,2,6,13,15,18,19)]

prot2_filtall_combo <- prot2_filtall %>% 
  group_by(KO_category2, newname1, treatment) %>% 
  summarize(tpm_sum=sum(TPM)) %>% 
  ungroup()

unique_genes <- unique(prot2_filtall_combo$newname1)
commB_summary <- aggregate(tpm_sum ~ newname1, data = prot2_filtall_combo[prot2_filtall_combo$treatment == "CommB", ], FUN = sum)  
sorted_genes <- commB_summary[order(commB_summary$tpm_sum, decreasing = FALSE), "newname1"]
genes_not_in_CommB <- unique_genes[!(unique_genes %in% sorted_genes)]
sorted_genes_all <- c(sorted_genes, genes_not_in_CommB)

prot2_filtall_combo$newname1 <- factor(prot2_filtall_combo$newname1, levels = sorted_genes_all)

prot2_filtall_noma <- filter(prot2_filtall_combo, ! KO_category2 == "Methionyl Aminopeptidase")#197
prot2_filtall_noma <- filter(prot2_filtall_noma, ! KO_category2 == "Other")#117

p2<- ggplot(prot2_filtall_noma,aes(x = treatment, y = newname1, fill = log2(tpm_sum+0.5))) +
  geom_tile() + 
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),  # Remove gridlines
    panel.background = element_rect(fill = "white", color = NA)) +
  facet_col(~KO_category2, scales = "free_y", space = "free")+theme(axis.text.y = element_blank())+theme(legend.position = "none")

# plot lfc
prot2_filt1 <- subset(tab1_filtered2, tab1_filtered2$gene_name %in% filtered_prot$gene_name)
prot2_filt2 <- subset(tab2_filtered2, tab2_filtered2$gene_name %in% filtered_prot$gene_name)
prot2_filt3 <- subset(tab3_filtered2, tab3_filtered2$gene_name %in% filtered_prot$gene_name)

prot2_filtall <- rbind(prot2_filt1, prot2_filt2, prot2_filt3)
prot2_filtall_def <- left_join(prot2_filtall, filtered_prot, by="gene_name")

prot2_filtall_def$newname1 <- paste(prot2_filtall_def$gene_name, " - ", prot2_filtall_def$KO_definition)

complete_data <- expand.grid(newname1 = unique(prot2_filtall_def$newname1),
                             treatment_comp = unique(prot2_filtall_def$treatment_comp))

# Merge complete dataset with original data to fill in missing logFC values
prot2_complete <- merge(complete_data, prot2_filtall_def, by = c("newname1", "treatment_comp"), all.x = TRUE)

# Replace missing logFC values with 0
prot2_complete$logFC[is.na(prot2_complete$logFC)] <- 0
prot2_complete$treatment_comp <- factor(prot2_complete$treatment_comp, levels=c("CommA", "CommB", "CommC"))
prot2_complete$newname1 <- factor(prot2_complete$newname1, levels = sorted_genes_all)

prot3_complete <- prot2_complete
prot3_complete <- subset(prot3_complete, KO_category == "Aminopeptidase")

# Apply the function to create the new column
prot3_complete <- mutate(prot3_complete, KO_category2 = sapply(KO_definition, find_category2))

prot3_noma <- filter(prot3_complete, ! KO_category2 == "Methionyl Aminopeptidase")#197
prot3_noma <- filter(prot3_noma, ! KO_category2 == "Other")#117

  
p1<- prot3_noma %>%   
  filter(!is.na(score))  %>% 
  ggplot( aes(x = logFC, y = newname1 , color = treatment_comp)) +
  geom_point(position = position_dodge(width = .5), size=4, alpha=.7) +  
  geom_vline(xintercept = 0, linetype = "dashed") +xlim(-5, 5)+
  scale_color_manual(values=c("#014b7a","orange", "#44b7c2")) +
  theme_minimal() +
  labs(x = "Log Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_col(~KO_category2, scales = "free_y", space = "free")+theme(axis.text.y = element_blank())+theme(legend.position = "none")


ggarrange(c2, c1, p2, p1, ncol=4, 
          widths = c(1, 2, 1, 2))


#### BACTERIAL PRODUCED PLANT GROWTH HORMONES ####
#IAA = "indol|tryptophan|phenylpyruvate|tyramine oxidase"
#pectinesterase:cell wall modification and growth
#cytokinins: zeatin | isopentenyl|dihydrozeatin
#ACC=non-proteinogenic amino acid ACC is the precursor and means of long-distance transport of ethylene, a plant hormone associated with growth arrest
#One of the most important genes associated with the increase in plant biomass and stress resistance is acdS, which encodes a 1-aminocyclopropane-1-carboxylate- or ACC-deaminase

#filter whole dataset by all important hormones
filtered_combi_horm <- separated_df %>%
  filter(str_detect(KO_definition, regex("aminocyclopropane|pectinesterase|zeatin|isopentenyl|dihydrozeatin|indol|tryptophan|phenylpyruvate|tyramine oxidase", ignore_case = TRUE)))

filtered_combi_horm  %>% 
  mutate(log2_TPM = log2(TPM+0.5)) %>%
  ggplot(aes(x = KO_definition, y = log2_TPM, group_by=treatment,fill = treatment)) +
  geom_boxplot(outlier.shape=NA, size=1, alpha=.7) + geom_jitter(aes(color=treatment),position = position_jitterdodge(0.2),cex=1) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_x_discrete(labels = function(KO_definition) str_wrap(KO_definition, width = 10))+ 
  scale_color_manual(values=col_treat)+facet_wrap(~week, ncol=1)+
  scale_fill_manual(values=col_treat)+ylab("log2(TPM+.5)")


# Plotting by LFC
filt_horm <- kegg.f2 %>%
  filter(str_detect(KO_definition, regex("aminocyclopropane|pectinesterase|indol|auxin|cytokinin|gibberellin|salicylic|zeatin|	
acetoin", ignore_case = TRUE)))

horm_filt1 <- subset(tab1_filtered2, tab1_filtered2$gene_name %in% filt_horm$gene_name)
horm_filt2 <- subset(tab2_filtered2, tab2_filtered2$gene_name %in% filt_horm$gene_name)
horm_filt3 <- subset(tab3_filtered2, tab3_filtered2$gene_name %in% filt_horm$gene_name)

horm_filtall <- rbind(horm_filt1, horm_filt2, horm_filt3)

horm_filtall$treatment_comp <- as.factor(horm_filtall$treatment_comp)
horm_filtall$treatment_comp <- factor(horm_filtall$treatment_comp, levels=c("CommC", "CommB", "CommA"))

horm_filtall_def <- left_join(horm_filtall, filt_horm, by="gene_name")
horm_filtall_def$newname1 <- paste(horm_filtall_def$gene_name, " - ", horm_filtall_def$KO_definition)
unique_genes <- unique(horm_filtall_def$newname1)
commB_summary <- aggregate(logFC ~ newname1, data = horm_filtall_def[horm_filtall_def$treatment_comp == "CommB", ], FUN = mean)  
sorted_genes <- commB_summary[order(commB_summary$logFC, decreasing = TRUE), "newname1"]
genes_not_in_CommB <- unique_genes[!(unique_genes %in% sorted_genes)]
sorted_genes_all <- c(sorted_genes, genes_not_in_CommB)
horm_filtall_def$newname1 <- factor(horm_filtall_def$newname1, levels = sorted_genes_all)
ggplot(horm_filtall_def, aes(x = newname1, y = treatment_comp, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0) +  # Adjust color scale as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### CARBOHYDRATE METABOLISM ####
separated_dfsorted  <- filter(separated_df, L2 %in% c("09101 Carbohydrate metabolism")) 

agg_df <- separated_dfsorted %>%
  group_by(treatment, KO_definition, KO) %>%
  summarize(total_TPM = sum(TPM, na.rm = TRUE))

agg_df  <-  agg_df %>% mutate(log2_TPM = log2(total_TPM+0.5))
commB_summary <- aggregate(log2_TPM ~ KO_definition, data = agg_df[agg_df$treatment == "CommB", ], FUN = mean)  
sorted_genes <- commB_summary[order(commB_summary$log2_TPM, decreasing = FALSE), "KO_definition"]
agg_df$KO_definition <- factor(agg_df$KO_definition, levels = sorted_genes)

detach("package:plyr", unload=TRUE)

agg_dffilt <- agg_df %>%
  group_by(KO_definition) %>%
  filter(
    # Check if CommB is highest compared to CommA or CommC
    any(log2_TPM[treatment == "CommB"] > log2_TPM[treatment == "CommA"]) |
      any(log2_TPM[treatment == "CommB"] > log2_TPM[treatment == "CommC"]),
    # Check if CommB is not lower than both CommA and CommC
    !(
      any(log2_TPM[treatment == "CommB"] < log2_TPM[treatment == "CommA"]) &
        any(log2_TPM[treatment == "CommB"] < log2_TPM[treatment == "CommC"])
    )
  )


agg_dffilt1 <- agg_dffilt %>% group_by(KO_definition) %>% 
  arrange(desc(log2_TPM),treatment="CommB") %>% 
  group_by(treatment) %>% slice(1:30)





agg_dffilt2 <- agg_df %>%
  group_by(KO_definition) %>%
  filter(
    # Condition 1: Keep rows where CommB is higher than both CommA and CommC
    all(log2_TPM[treatment == "CommB"] >= log2_TPM[treatment == "CommA"]) ||
      all(log2_TPM[treatment == "CommB"] >= log2_TPM[treatment == "CommC"]),
    
    # Condition 2: Remove rows where CommB is only marginally higher than CommA or CommC
    all(log2_TPM[treatment == "CommB"] - log2_TPM[treatment == "CommA"] >= 4) ||
      all(log2_TPM[treatment == "CommB"] - log2_TPM[treatment == "CommC"] >= 4)
  )


ecoc<- ggplot(agg_dffilt2,aes(x = treatment, y = KO_definition, fill = log2_TPM)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0, limits = c(NA, 12)) +  # Adjust color scale as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "09101 Carbohydrate metabolism")




#### AMINO ACID METABOLISM ####
separated_dfsorted  <- filter(separated_df, L2 %in% c("09105 Amino acid metabolism")) 

agg_df <- separated_dfsorted %>%
  group_by(treatment, KO_definition) %>%
  summarize(total_TPM = sum(TPM, na.rm = TRUE))

agg_df  <-  agg_df %>% mutate(log2_TPM = log2(total_TPM+0.5))
commB_summary <- aggregate(log2_TPM ~ KO_definition, data = agg_df[agg_df$treatment == "CommB", ], FUN = mean)  
sorted_genes <- commB_summary[order(commB_summary$log2_TPM, decreasing = FALSE), "KO_definition"]
agg_df$KO_definition <- factor(agg_df$KO_definition, levels = sorted_genes)

agg_dffilt <- agg_df %>%
  group_by(KO_definition) %>%
  filter(
    # Check if CommB is highest compared to CommA or CommC
    any(log2_TPM[treatment == "CommB"] > log2_TPM[treatment == "CommA"]) |
      any(log2_TPM[treatment == "CommB"] > log2_TPM[treatment == "CommC"]),
    # Check if CommB is not lower than both CommA and CommC
    !(
      any(log2_TPM[treatment == "CommB"] < log2_TPM[treatment == "CommA"]) &
        any(log2_TPM[treatment == "CommB"] < log2_TPM[treatment == "CommC"])
    )
  )


agg_dffilt2 <- agg_df %>%
  group_by(KO_definition) %>%
  filter(
    # Condition 1: Keep rows where CommB is higher than both CommA and CommC
    all(log2_TPM[treatment == "CommB"] >= log2_TPM[treatment == "CommA"]) ||
      all(log2_TPM[treatment == "CommB"] >= log2_TPM[treatment == "CommC"]),
    
    # Condition 2: Remove rows where CommB is only marginally higher than CommA or CommC
    all(log2_TPM[treatment == "CommB"] - log2_TPM[treatment == "CommA"] >=4) ||
      all(log2_TPM[treatment == "CommB"] - log2_TPM[treatment == "CommC"] >= 4)
  )

ecoa<- agg_dffilt2  %>% 
  ggplot(aes(x = treatment, y = KO_definition, fill = log2_TPM)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3498DB", mid = "white", high = "#E74C3C", midpoint = 0, limits = c(NA, 12)) +  # Adjust color scale as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("09105 Amino acid metabolism")


ggarrange(ecoc, ecoa, nrow=1, widths=c(1,1.15))

## EcoPlate LFC
filt_eco <- KEGG_KO2 %>%
  filter(str_detect(KO_definition, regex("methyltransferases|cellobiose|GlcNAc|acetylglucosamine|glucosidase", ignore_case = TRUE)))

filt_eco1 <- subset(tab1_filtered2, tab1_filtered2$gene_name %in% filt_eco$gene_name)
filt_eco2 <- subset(tab2_filtered2, tab2_filtered2$gene_name %in% filt_eco$gene_name)
filt_eco3 <- subset(tab3_filtered2, tab3_filtered2$gene_name %in% filt_eco$gene_name)

filt_ecoall <- rbind(filt_eco1, filt_eco2, filt_eco3)

filt_ecoall$treatment_comp <- as.factor(filt_ecoall$treatment_comp)
filt_ecoall$treatment_comp <- factor(filt_ecoall$treatment_comp, levels=c("CommC", "CommB", "CommA"))
eco_filtall_def <- left_join(filt_ecoall, KEGG_KO2, by="gene_name")


eco_filtall_def$newname1 <- paste(eco_filtall_def$gene_name, " - ", eco_filtall_def$KO_definition)
commB_summary <- aggregate(logFC ~ newname1, data = eco_filtall_def[eco_filtall_def$treatment_comp == "CommB", ], FUN = mean)  
sorted_genes <- commB_summary[order(commB_summary$logFC, decreasing = TRUE), "newname1"]
genes_not_in_CommB <- unique_genes[!(unique_genes %in% sorted_genes)]
sorted_genes_all <- c(sorted_genes, genes_not_in_CommB)
eco_filtall_def$newname1 <- factor(eco_filtall_def$newname1, levels = sorted_genes_all)
a <- eco_filtall_def %>% filter(L1 == "09100 Metabolism")

ggplot(a, aes(x = newname1, y = treatment_comp, fill = logFC)) +
  geom_tile() +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick") +  # Adjust color scale as needed
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6))


#### NMDS FOR FUNCTIONS ####
dim(tpm_fun_subset)
#162623
tpm2 <- as.data.frame(tpm_fun_subset)

#transform
tpm3t <- t(tpm2)
rownames(meta_rna2) <- meta_rna2$sample_id
#match up with metadata
rownames(tpm3t) == rownames(meta_rna2)
meta_rna2$treatment <- factor(meta_rna2$treatment, levels = c("CommA", "CommB", "CommC"))

# Arrange the dataframe by 'treatment' and 'week'
meta_rna2 <- meta_rna2 %>%
  arrange(treatment, week)
tpm3t <- as.data.frame(tpm3t)
matched_cols <- rownames(tpm3t)[match(rownames(meta_rna2), rownames(tpm3t))]
tpm3t <- tpm3t[matched_cols,]
rownames(tpm3t) == rownames(meta_rna2)

#run NMDS on genes
set.seed(123)
nmds1 <- metaMDS(tpm3t, k=2, trymax=1000)
#Dimensions: 2 
#Stress:     0.09904656   
#Stress type 1, weak tie

all.scores <- as.data.frame(scores(nmds1, "sites"))
#add metadata to the scores info
all.scores <- cbind(all.scores, meta_rna2)

#add hulls
all.scores <- all.scores %>% 
  unite(hull.id, treatment, remove = FALSE)

hullsall <- all.scores %>%
  group_by(hull.id) %>%
  slice(chull(NMDS1, NMDS2))

ggplot(all.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = treatment, shape = week)) +
  theme(axis.text.y = element_text(colour = "black", size = 16, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 16), 
        legend.text = element_text(size = 16, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, colour = "black"), 
        legend.title = element_text(size = 16, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "treatment", y = "NMDS2", shape = "week")  + 
  scale_colour_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))  + aes(fill = factor(hull.id)) + geom_polygon(data = hullsall, alpha = 0.2) + guides(fill=FALSE)+ 
  scale_fill_manual(values = c("#014b7a", "#ffaf49", "#44b7c2"))



# set number of permutations
perm <- how(nperm = 999)
adonis <- adonis2(tpm3t ~ treatment * week, data=all.scores, permutations = perm, method = "bray")
adonis


#### VENN DIAGRAM ####
treatmentA_genes <- tab_filtered2$gene_name[tab_filtered2$treatment_comp == "CommA"]
treatmentB_genes <- tab_filtered2$gene_name[tab_filtered2$treatment_comp == "CommB"]
treatmentC_genes <- tab_filtered2$gene_name[tab_filtered2$treatment_comp == "CommC"]

gene_lists <- list(
  CommA = treatmentA_genes,
  CommB = treatmentB_genes,
  CommC = treatmentC_genes
)

venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("CommA", "CommB", "CommC"),
  filename = NULL,  # Set to NULL to plot directly in R
  output = FALSE
)
dev.off()
grid.draw(venn.plot)


#### TOP 20 HEATMAP ####
summary_table <- tab_filtered2 %>%
  mutate(status = case_when(
    logFC > 0 & adj.P.Val < 0.05 ~ "up",
    logFC < 0 & adj.P.Val < 0.05 ~ "down",
    TRUE ~ "not.significant"
  )) %>%
  group_by(treatment_comp, status) %>%
  summarize(count = n()) %>%
  filter(status != "not.significant")

summary_table$status <- factor(summary_table$status, levels=c("up", "down"))

# Plot the bar chart
ggplot(summary_table, aes(x = treatment_comp, y = count, fill = status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Differential Gene Expression Summary", x = "Treatment", y = "Gene Count") +
  scale_fill_manual(values = c("up" = "#E74C3C", "down" = "#3498DB")) +
  theme_minimal()

tab_filtered2ma <- top_all %>%
  mutate(status = case_when(
    adj.P.Val < 0.05 & logFC > 0 ~ "upregulated",
    adj.P.Val < 0.05 & logFC < 0 ~ "downregulated",
    TRUE ~ "not.significant"
  ))

tab_filtered2ma$status <- factor(tab_filtered2ma$status, levels = c("not.significant", "upregulated", "downregulated"))

# Create MA plots for each treatment
for (treatment in unique(tab_filtered2ma$treatment_comp)) {
  df_treatment <- tab_filtered2ma %>% filter(treatment_comp == treatment)
  
  p <- ggplot() +
    geom_point(data = df_treatment %>% filter(status == "not.significant"), aes(x = AveExpr, y = logFC), shape=16, color = "gray", alpha = 0.6) +
    geom_point(data = df_treatment %>% filter(status == "upregulated"), aes(x = AveExpr, y = logFC), shape=16, color = "#E74C3C", alpha = 0.6) +
    geom_point(data = df_treatment %>% filter(status == "downregulated"), aes(x = AveExpr, y = logFC), shape=16, color = "#3498DB", alpha = 0.6) +
    labs(title = paste("MA Plot -", treatment), x = "Average Expression", y = "Log Fold Change") +
    theme_minimal() +
    theme(legend.position = "right")+ylim(c(-11, 11))
  
  print(p)  # Print the plot to display it
}


#### separate heat map for each treatment
kegg_filt <- KEGG_KO[,c(1,6)]
kegg_dist <- kegg_filt %>%
  distinct()

common_theme <- theme(
  legend.position = "none",
  text = element_text(size = 10)  # Adjust font size here
)

# Filter by treatment and select the top 20 highest LFC genes for CommA
ave_fun <- tab_filtered2
ave_fun_kegg <- ave_fun %>%
  left_join(kegg_dist, by = "gene_name")

ave_fun_kegg2 <- ave_fun_kegg %>%
  group_by(treatment_comp, KO_definition) %>%
  summarize(ave_lfc = mean(logFC, na.rm = TRUE)) %>%
  ungroup()
detach("package:plyr", unload=TRUE)
top_ave_a <- ave_fun_kegg2 %>%
  filter(treatment_comp == "CommA") %>%
  arrange(desc(ave_lfc)) %>%
  slice_head(n = 15)

top_ave_a <- top_ave_a %>%
  mutate(KO_definition = make.unique(as.character(KO_definition))) %>%
  mutate(KO_definition = factor(KO_definition, levels = KO_definition[order(ave_lfc)]))

ko_a <- top_ave_a$KO_definition
filtered_dfa <- kegg_filt[kegg_filt$KO_definition %in% ko_a, ]

library(viridis)

d<-   ggplot(top_ave_a, aes(x = treatment_comp, y = KO_definition, fill = ave_lfc)) +
geom_tile() +
  scale_fill_viridis(option = "C", limits = c(3, 8)) + # "C" option gives a nice gradient for subtle differences
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ common_theme

###COMMB
top_ave_b <- ave_fun_kegg2 %>%
  filter(treatment_comp == "CommB") %>%
  arrange(desc(ave_lfc)) %>%
  slice_head(n =15)

top_ave_b <- top_ave_b %>%
  mutate(KO_definition = make.unique(as.character(KO_definition))) %>%
  mutate(KO_definition = factor(KO_definition, levels = KO_definition[order(ave_lfc)]))

ko_b <- top_ave_b$KO_definition
filtered_dfb <- kegg_filt[kegg_filt$KO_definition %in% ko_b, ]


e<- ggplot(top_ave_b, aes(x = treatment_comp, y = KO_definition, fill = ave_lfc)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(3, 8)) + # "C" option gives a nice gradient for subtle differences
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ common_theme


####COMMC
top_ave_c <- ave_fun_kegg2 %>%
  filter(treatment_comp == "CommC") %>%
  arrange(desc(ave_lfc)) %>%
  slice_head(n = 15)

top_ave_c <- top_ave_c %>%
  mutate(KO_definition = make.unique(as.character(KO_definition))) %>%
  mutate(KO_definition = factor(KO_definition, levels = KO_definition[order(ave_lfc)]))


ko_c <- top_ave_c$KO_definition
filtered_dfc <- kegg_filt[kegg_filt$KO_definition %in% ko_c, ]

f<- ggplot(top_ave_c, aes(x = treatment_comp, y = KO_definition, fill = ave_lfc)) +
  geom_tile() +
  scale_fill_viridis(option = "C", limits = c(3, 8)) + # "C" option gives a nice gradient for subtle differences
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ common_theme




ggarrange(d, e, f, 
          ncol = 3, nrow = 1,widths = c(.5, .5, .5))


#### MEAN ABUNDANCE PLOTS ####
### CommA
ta <- topA2 %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))

#make rowname a column called gene_name
ta$gene_name <- rownames(ta)
ta$delabel <- NA
ta$delabel[ta$status != "not.signif"] <- rownames(ta)[ta$status != "not.signif"]
ta$outliers <- NA



outlier_conditions <- ta$gene_name %in% filtered_dfa$gene_name

# Update the outliers column using the conditions
ta$outliers[outlier_conditions & ta$status != "not.signif"] <- rownames(ta)[outlier_conditions & ta$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
ta$gene_name <- as.factor(ta$gene_name)
tpm_fun_meta$gene_name <- as.factor(tpm_fun_meta$gene_name)

ta_modified_KO <- left_join(ta, tpm_fun_meta, by = "gene_name")
na_rows <- is.na(ta_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
ta_modified_KO[na_rows, 12:16] <- NA
ta_modified_KO <- ta_modified_KO %>% arrange(status)

#FIG7
ggplot(ta_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(shape=20,size=2, alpha = .5, stroke=NA) +
  scale_color_manual(values=c("grey","#E74C3C",  "#3498DB")) +
  theme_classic()+ 
  geom_text_repel(aes(label = KO_definition), size = 2, color="black", box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)+
  theme(legend.position = "none")

### CommB
tb <- topB2 %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
tb$gene_name <- rownames(tb)
tb$delabel <- NA
tb$delabel[tb$status != "not.signif"] <- rownames(tb)[tb$status != "not.signif"]
tb$outliers <- NA

outlier_conditions <- tb$gene_name %in% filtered_dfb$gene_name

# Update the outliers column using the conditions
tb$outliers[outlier_conditions & tb$status != "not.signif"] <- rownames(tb)[outlier_conditions & tb$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
tb_modified_KO <- left_join(tb, tpm_fun_meta, by = "gene_name")
na_rows <- is.na(tb_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
tb_modified_KO[na_rows, 12:16] <- NA
tb_modified_KO <- tb_modified_KO %>% arrange(status)

#FIG7
ggplot(tb_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(shape=20,size=2, alpha = .5, stroke=NA) +
  scale_color_manual(values=c("grey","#E74C3C",  "#3498DB")) +
  theme_classic()+ 
  geom_text_repel(aes(label = KO_definition), size = 2, color="black", box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)+
  theme(legend.position = "none")

### CommC
tc <- topC2 %>% 
  mutate(status=factor(case_when(logFC>0 & adj.P.Val<0.05 ~ "up",
                                 logFC<0 & adj.P.Val<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("not.signif","up","down")))
#make rowname a column called gene_name
tc$gene_name <- rownames(tc)
tc$delabel <- NA
tc$delabel[tc$status != "not.signif"] <- rownames(tc)[tc$status != "not.signif"]
tc$outliers <- NA

outlier_conditions <- tc$gene_name %in% filtered_dfc$gene_name

# Update the outliers column using the conditions
tc$outliers[outlier_conditions & tc$status != "not.signif"] <- rownames(tc)[outlier_conditions & tc$status != "not.signif"]

#combine the kegg.f so we can see what the function of the outliers are
tc_modified_KO <- left_join(tc, tpm_fun_meta, by = "gene_name")
na_rows <- is.na(tc_modified_KO$outliers)

# Set columns 12 to 20 to NA for the identified rows
tc_modified_KO[na_rows, 12:16] <- NA
tc_modified_KO <- tc_modified_KO %>% arrange(status)

#FIG7
ggplot(tc_modified_KO, aes(x=AveExpr, y=logFC, color=status)) +
  geom_point(shape=20,size=2, alpha = .5, stroke=NA) +
  scale_color_manual(values=c("grey","#E74C3C",  "#3498DB")) +
  theme_classic()+ 
  geom_text_repel(aes(label = KO_definition), size = 2,color="black", box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf)+
  theme(legend.position = "none")


