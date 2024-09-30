# Packages
library(tidyverse)
library(ComplexUpset)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(mixOmics)
library(vegan)
library(caret)
library(limma)
library(patchwork)
library(rstatix)

########################################################
# FEATURE TABLE BLANK SUBTRACTION - POSITIVE IONIZATION#
########################################################

metadata <- read.csv("Metadata/3DMM_metadata.csv")
feature_pos <- read.csv(file = "Feature_tables/3DMM_pos_quant.csv")
colnames(feature_pos) <- gsub(".mzML.Peak.area", "", colnames(feature_pos))
metadata_filtered <- metadata[metadata$Well %in% colnames(feature_pos),]

# For each plate, check if sample peak area is five times higher than the average of blanks and impute as zero if not.
plate <- c("P1_","P2_","P3_","P4_","P5_","P6_","P7_","P8_","P9_","P10_")
feature_pos_subtract <- feature_pos[,c("row.ID", "row.m.z", "row.retention.time")]
for (i in 1:length(plate)){
  metadata_plate <- metadata_filtered[grepl(plate[i], metadata_filtered$Well),]  # map for specific plate
  swabblk_plate <- metadata_plate$Well[metadata_plate$Barcode == "swab.blank"]  # name vector for swab blk in this plate
  sample_plate <- metadata_plate$Well[!metadata_plate$Barcode %in% c("swab.blank","solvent.blank")]  # name vector for sample in this plate
  swabblk_avr <- rowMeans(feature_pos[,swabblk_plate]) # average of swab blk
  feature_sample_plate <- feature_pos[,sample_plate] 
  t <- (feature_sample_plate > 5*swabblk_avr)  # if higher than 5*blks
  feature_sample_plate_subtracted <- t * feature_sample_plate - swabblk_avr  #subtract blk
  feature_sample_plate_subtracted <- (feature_sample_plate_subtracted > 0) * feature_sample_plate_subtracted
  feature_pos_subtract <- cbind(feature_pos_subtract, feature_sample_plate_subtracted)
}

# remove feature with less than two detections in all modules
feature_pos_filtered <- feature_pos_subtract %>% column_to_rownames("row.ID") %>% dplyr::select(contains("_")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID")
feature_metadata <- merge(feature_pos_filtered, metadata_filtered[,c("Well", "Module")], 
                          by.x = "SampleID", by.y = "Well", all.x = T, all.y = F)
feature_module_count <- feature_metadata %>% column_to_rownames("SampleID") %>% 
  aggregate(list(feature_metadata$Module), function(c)sum(c!=0)) %>% 
  dplyr::select(!contains("Module")) %>% column_to_rownames("Group.1")
feature_keep <- colnames(feature_module_count)[apply(feature_module_count, 2, function(x) any(x > 1))]

# save feature module count table for future processing
feature_module_count_filtered <- feature_module_count %>% 
  dplyr::select(all_of(feature_keep)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("feature.ID")
write.csv(feature_module_count_filtered, "Feature_tables/3DMM_pos_feature_module_count.csv", row.names = F)

# save blank removed feature table for FBMN job on GNPS
feature_pos_duplicate_filtered <-  feature_pos_subtract %>% dplyr::filter(feature_pos_subtract$row.ID %in% feature_keep)
# modify the feature table format to be compatible with GNPS FBMN workflow
colnames(feature_pos_duplicate_filtered)[1:3] <- c("row ID", "row m/z", "row retention time")
colnames(feature_pos_duplicate_filtered)[4:ncol(feature_pos_duplicate_filtered)] <- 
  paste(colnames(feature_pos_duplicate_filtered)[4:ncol(feature_pos_duplicate_filtered)], ".mzML Peak area", sep = "")
# export feature table
write.csv(feature_pos_duplicate_filtered, "Feature_tables/3DMM_pos_quant_blk_duplicate_filtered.csv", row.names = F)

#################################################################
# FILTER MGF FILE TO REMOVE BLANK FEATURES - POSITIVE IONIZATION#
#################################################################

# define a function to only keep features in the feature_keep for a mgf
clean_mgf_featureID <- function(input_file, feature_to_keep){
  blocks_to_keep <- list() # Initialize an empty list to store the spectra to keep
  current_block <- "" # Initialize an empty string to store the current block
  for (line in input_file) {
    if(current_block == ""){
      current_block <- line # add the line to the current block and avoid an empty line before the spectra
    }
    else{
      current_block <- paste(current_block, line, sep = "\n")
    }
    if (startsWith(line, "END IONS")) {
      feature_id_line <- grep("FEATURE_ID", unlist(strsplit(current_block, "\n")), value = TRUE)
      feature_id <- as.numeric(sub("FEATURE_ID=", "", feature_id_line)) # Extract FEATURE_ID from the current block
      if (feature_id %in% feature_to_keep) {
        blocks_to_keep <- c(blocks_to_keep, current_block) # Append the current block to blocks_to_keep if FEATURE_ID is in keep list
      }
      current_block <- ""  # Reset current_block for the next block
    }
  }
  cleaned_content <- paste(blocks_to_keep, collapse = "\n\n") # Combine the blocks to keep into a single string with an empty line between each block
  return(cleaned_content)
}

# define a function to remove doulbe MS1 scan for some features from the SIRIUS mgf file
spectype_remove <- function(input_file){
  blocks_to_keep <- list() # Initialize an empty list to store the spectra to keep
  current_block <- "" # Initialize an empty string to store the current block
  for (line in input_file) {
    if(current_block == ""){
      current_block <- line # add the line to the current block and avoid an empty line before the spectra
    }
    else{
      current_block <- paste(current_block, line, sep = "\n")
    }
    if (startsWith(line, "END IONS")) {
      if (sum(grepl("SPECTYPE", unlist(strsplit(current_block, "\n")))) == 0) {
        blocks_to_keep <- c(blocks_to_keep, current_block) # Append the current block to blocks_to_keep if FEATURE_ID is in keep list
      }
      current_block <- ""  # Reset current_block for the next block
    }
  }
  cleaned_content <- paste(blocks_to_keep, collapse = "\n\n") # Combine the blocks to keep into a single string with an empty line between each block
  return(cleaned_content)
}

# modify mgf file to remove blank features
mgf_raw <- readLines("MGF/3DMM_pos_FBMN.mgf") 
mgf_modified <- clean_mgf_featureID(mgf_raw, feature_keep)
writeLines(mgf_modified, "MGF/3DMM_pos_FBMN_blk_rm.mgf") # mgf file for FBMN job on GNPS and fastMASST analysis
# sirius mgf
mgf_sirius_raw <- readLines("MGF/3DMM_pos_SIRIUS.mgf")
mgf_sirius_modified <- clean_mgf_featureID(mgf_sirius_raw, feature_keep)
writeLines(mgf_sirius_modified, "MGF/3DMM_pos_SIRIUS_blk_rm.mgf") 
mgf_sirius_blk_rm <- readLines("MGF/3DMM_pos_SIRIUS_blk_rm.mgf")
mgf_sirius_spectype_removed <- spectype_remove(mgf_sirius_blk_rm)
writeLines(mgf_sirius_spectype_removed, "MGF/3DMM_pos_SIRIUS_blk_spectype_rm.mgf") # mgf file for SIRIUS

######################################################
# UPSET PLOT TO VISUALIZE FEATURE MODULE DISTRIBUTION# 
######################################################

# SIRIUS, Qemistree, and qiime2 for alpha-diversity and beta-diversity were run by Dr. Helena Russo with code in Qemistree_Space_Station.ipynb
# Output result files are stored in folder "Qiime_artifact"

# UpSet plot for feature distribution in different modules, stacked bar plot colored by structural class prediction
feature_module_count <- read.csv("Feature_tables/3DMM_pos_feature_module_count.csv")
feature_class <- read.csv("Qiime_artifact/3DMM_pos_feature_structure_class.csv")
feature_combined <- merge(feature_class, feature_module_count, by.x = "X.featureID", by.y = "feature.ID", all = T)
structure_class_count <- as.data.frame(table(feature_combined$ClassyFire.superclass))
class_ignore <- structure_class_count$Var1[structure_class_count$Freq < 20]
feature_combined$ClassyFire.superclass_cleaned <- feature_combined$ClassyFire.superclass
feature_combined$ClassyFire.superclass_cleaned[is.na(feature_combined$ClassyFire.superclass_cleaned) | 
                                                 feature_combined$ClassyFire.superclass_cleaned %in% class_ignore] <- "Not classified"
feature_combined$ClassyFire.superclass_cleaned <- factor(feature_combined$ClassyFire.superclass_cleaned, levels = c(
  "Benzenoids", "Lignans, neolignans and related compounds", "Lipids and lipid-like molecules",
  "Nucleosides, nucleotides, and analogues", "Organic acids and derivatives", "Organic nitrogen compounds",
  "Organic oxygen compounds", "Organohalogen compounds", "Organoheterocyclic compounds", 
  "Phenylpropanoids and polyketides", "Not classified" 
))
module <- c("Airlock", "Columbus", "JLP", "JPM", "Node.1", "Node.2", "Node.3", "PMM", "US.Lab")
feature_combined[module] <- (feature_combined[module] > 0) # convert module occurrence count to True/False

#plot module unique feature
upset(feature_combined, module, 
      base_annotations=list(
        'Intersection'=intersection_size(
          counts=TRUE, bar_number_threshold = 500, 
          mapping=aes(fill=ClassyFire.superclass_cleaned, color = "black"))
        + scale_fill_manual(values=c("Lipids and lipid-like molecules"="#CC0066",
                                     "Organic oxygen compounds"="#66CC99",
                                     "Organic acids and derivatives"="#FFFF99",
                                     "Benzenoids"="#666600",
                                     "Nucleosides, nucleotides, and analogues"="#CCCCCC",
                                     "Organoheterocyclic compounds"="#FF9966",
                                     "Organic nitrogen compounds"="#996699",
                                     "Organohalogen compounds"="#6699CC",
                                     "Phenylpropanoids and polyketides"="#FF99CC",
                                     "Lignans, neolignans and related compounds"="#0099FF",
                                     "Not classified"="black"))
      ), 
      max_degree=2, 
      height_ratio=0.3,
      sort_intersections_by=c('degree', 'cardinality'),
)

#plot common feature
upset(feature_combined, module, 
      base_annotations=list(
        'Intersection'=intersection_size(
          counts=TRUE, bar_number_threshold = 500, 
          mapping=aes(fill=ClassyFire.superclass_cleaned, color = "black"))
        + scale_fill_manual(values=c("Lipids and lipid-like molecules"="#CC0066",
                                     "Organic oxygen compounds"="#66CC99",
                                     "Organic acids and derivatives"="#FFFF99",
                                     "Benzenoids"="#666600",
                                     "Nucleosides, nucleotides, and analogues"="#CCCCCC",
                                     "Organoheterocyclic compounds"="#FF9966",
                                     "Organic nitrogen compounds"="#996699",
                                     "Organohalogen compounds"="#6699CC",
                                     "Phenylpropanoids and polyketides"="#FF99CC",
                                     "Lignans, neolignans and related compounds"="#0099FF",
                                     "Not classified"="black"))
      ), 
      min_degree=7, min_size = 50,
      height_ratio=0.3,
      sort_intersections_by=c('cardinality'),
)

####################################
# HIERARCY-INFORMED ALPHA-DIVERSITY# 
####################################
# Faith's phylogenetic diversity metric was calculated in qiime2  by Dr. Helena Russo
# code: Qemistree_Space_Station.ipynb. Result table in folder "Qiime_artifact"

mb_alpha <- read.csv("Qiime_artifact/3DMM_pos_alpha_diversity_faith_pd.csv")
metadata <- read.csv("Metadata/3DMM_metadata.csv")
mb_alpha$X.SampleID <- gsub(".mzML", "", mb_alpha$X.SampleID)
alpha_metadata <- merge(mb_alpha, metadata, by.x="X.SampleID", by.y="Well")
alpha_metadata <- alpha_metadata[alpha_metadata$X.SampleID != "P3_G5",]

# Boxplot to visualize module distribution
module_ordered <- with(alpha_metadata, reorder(Module, faith_pd, median, decreasing = F))
alpha_metadata$Module <- factor(alpha_metadata$Module, levels = levels(module_ordered))
ggplot(alpha_metadata, aes(x=Module, y=faith_pd, color="black", fill=Module)) + 
  geom_boxplot(outlier.shape=NA, color = "black", linewidth=1) + #avoid plotting outliers twice
  stat_boxplot(geom = "errorbar", width = 0.5, color="black", linewidth=1) + 
  geom_jitter(position=position_jitter(width=.1, height=0), size = 1.5, col = "black")+
  scale_fill_manual(values = c("#5DA09F","#0E8040", "#F27521","#F3EC1A","#3954A5","#92298D","#ED2024","#F39AC1","#8ED5E3")) +
  theme_classic()+
  labs(x = NULL, y = "Metabolite Faith's Phylo-Diversity") +
  theme(axis.ticks.length = unit(1, "mm"), 
        axis.ticks.x.bottom = element_line(linewidth=0.8),
        axis.ticks.y.left = element_line(linewidth=0.8),
        axis.line.y.left = element_line(size = 0.8),
        axis.line.x.bottom=element_line(linewidth=0.8),
        axis.text.x = element_text(size = 21, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 21, color = "black"),
        axis.title = element_text(size = 21),
        legend.position = "none"
  )
oneway.test(faith_pd ~ Plate, data = alpha_metadata, var.equal=T)

# Correlation with genomics alpha-diversity 
seq_alpha <- read.csv("Qiime_artifact/3DMM_genomics_alpha_diversity_table.csv")
alpha_metadata$Labels <- sub("-",".",alpha_metadata$Labels)
alpha_metadata$Labels <- sub("/",".",alpha_metadata$Labels)
rownames(alpha_metadata) <- NULL
alpha_metadata_filtered <- alpha_metadata %>% dplyr::select(Labels, faith_pd) %>% column_to_rownames("Labels")
colnames(alpha_metadata_filtered) <- "faith_pd_metabolomics"
seq_alpha_filtered <- seq_alpha %>% dplyr::filter(seq_alpha$SampleID %in% rownames(alpha_metadata_filtered)) %>% column_to_rownames("SampleID")
df_merged <- merge(seq_alpha_filtered, alpha_metadata_filtered, by = "row.names", all.x = T, all.y = F)
# Generate scatter plot with regression line
sp <- ggscatter(df_merged, x = "faith_pd_metabolomics", y = "faith_pd",
                color = "#0073C2FF", alpha = 0.6, size = 3.5,
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Metabolite Faith's Phylo-Diversity", ylab = "16s Faith's Phylo-Diversity") +
  theme(axis.text.y = element_text(size=21),
        axis.text.x = element_text(size=21),
        axis.title.y = element_text(size=21),
        axis.title.x = element_text(size=21),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(1, "mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(df_merged, "faith_pd_metabolomics", fill = "gray")
yplot <- ggdensity(df_merged, "faith_pd", fill = "gray")+rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
p1 <- insert_xaxis_grob(sp, xplot, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right") 
ggdraw(p2)

####################################################################
# CORRELATE DISINFECTANT INTENSITY WITH MICROBIOME ALPHA-DIVERSITY # 
####################################################################
# Disinfectant feature:
# Quaternary ammonium compounds: 16552, 18230, 18052, 20834; Transformation products: 12752, 13456, 11247, 13982, 14285
# Cocamidopropyl betaine derivatives: 12470, 15143

seq_alpha <- read.csv("Qiime_artifact/3DMM_genomics_alpha_diversity_table.csv")
chem_feature <- read.csv("Feature_tables/3DMM_pos_quant_blk_duplicate_filtered.csv")
metadata <- read.csv("Metadata/3DMM_metadata.csv")
colnames(chem_feature) <- sub(".mzML.Peak.area", "", colnames(chem_feature))

# Convert the sample name from vial location to sample label
chem_feature_filtered <- chem_feature %>%
  column_to_rownames("row.ID") %>%
  dplyr::select(contains("_")) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Sample.ID") %>%
  left_join(metadata, by = c("Sample.ID" = "Well")) %>%
  column_to_rownames("Labels") %>%
  dplyr::select(all_of(as.character(chem_feature$row.ID))) 
colnames(chem_feature_filtered) <- paste("F", colnames(chem_feature_filtered), sep = "")
rownames(chem_feature_filtered) <- sub("-",".",rownames(chem_feature_filtered))
rownames(chem_feature_filtered) <- sub("/",".",rownames(chem_feature_filtered))
seq_alpha_filtered <- seq_alpha %>% dplyr::filter(seq_alpha$SampleID %in% rownames(chem_feature_filtered)) %>% column_to_rownames("SampleID")
df_merged <- merge(seq_alpha_filtered, chem_feature_filtered, by = "row.names", all.x = T, all.y = F)

#group different disinfactant group
data_sum <- data.frame(faith_pd = df_merged$faith_pd, 
                       sum_QAC_betaine = (df_merged$F16552+df_merged$F18230+df_merged$F18052+df_merged$F20834
                                          +df_merged$F12752+df_merged$F13456+df_merged$F11247+df_merged$F13982+df_merged$F14285
                                          +df_merged$F12470+df_merged$F15143))
data_log <- cbind(data_sum, log_QAC_bt = log(data_sum$sum_QAC_betaine+1)) # Log-transform the metabolomics peak area
data_log_filtered <- data_log[data_log$log_QAC_bt >0,]

# Scatter plot with regression line
sp <- ggscatter(data_log_filtered, x = "log_QAC_bt", y = "faith_pd", 
                color = "brown", alpha = 0.6, size = 3.5,
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "Disinfectants Peak Area", ylab = "16S Faith's Phylo-Diversity") +
  theme(axis.text.y = element_text(size=21),
        axis.text.x = element_text(size=21),
        axis.title.y = element_text(size=21),
        axis.title.x = element_text(size=21),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(1, "mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(data_log_filtered, "log_QAC_bt", fill = "gray")
yplot <- ggdensity(data_log_filtered, "faith_pd", fill = "gray")+rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
p1 <- insert_xaxis_grob(sp, xplot, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right") 
ggdraw(p2)

#########################################################
# CORRELATE DISINFECTANT INTENSITY WITH 16S READ COUNTS # 
#########################################################
unrarefied <- read.csv("16s_data/katharo_deblur_gg2_unrarefied_genera_collapsed.csv", row.names = 1)
rownames(unrarefied) <- sub("151322.14542.", "", rownames(unrarefied))
rownames(unrarefied) <- sub("150297.14542.", "", rownames(unrarefied))

chem_feature <- read.csv("Feature_tables/3DMM_pos_quant_blk_duplicate_filtered.csv")
metadata <- read.csv("Metadata/3DMM_metadata.csv")
colnames(chem_feature) <- sub(".mzML.Peak.area", "", colnames(chem_feature))

# Convert the sample name from vial location to sample label
chem_feature_filtered <- chem_feature %>%
  column_to_rownames("row.ID") %>%
  dplyr::select(contains("_")) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("Sample.ID") %>%
  left_join(metadata, by = c("Sample.ID" = "Well")) %>%
  column_to_rownames("Labels") %>%
  dplyr::select(all_of(as.character(chem_feature$row.ID)))
colnames(chem_feature_filtered) <- paste("F", colnames(chem_feature_filtered), sep = "")
rownames(chem_feature_filtered) <- sub("-",".",rownames(chem_feature_filtered))
rownames(chem_feature_filtered) <- sub("/",".",rownames(chem_feature_filtered))
unrarefied_filtered <- unrarefied %>% dplyr::filter(rownames(unrarefied) %in% rownames(chem_feature_filtered))
unrarefied_filtered <- unrarefied_filtered %>% mutate(counts = log10(rowSums(.)))
df_merged <- merge(unrarefied_filtered %>% dplyr::select(counts), chem_feature_filtered, by = "row.names", all.x = T, all.y = F)

#group different disinfectant group
data_sum <- data.frame(counts = df_merged$counts, 
                       sum_QAC_betaine = (df_merged$F16552+df_merged$F18230+df_merged$F18052+df_merged$F20834
                                          +df_merged$F12752+df_merged$F13456+df_merged$F11247+df_merged$F13982+df_merged$F14285
                                          +df_merged$F12470+df_merged$F15143))
data_log <- cbind(data_sum, log_QAC_bt = log(data_sum$sum_QAC_betaine+1)) # Log-transform the metabolomics peak area
data_log_filtered <- data_log[data_log$log_QAC_bt >0,]

# Scatter plot with regression line
sp <- ggscatter(data_log_filtered, x = "log_QAC_bt", y = "counts", 
                color = "orange", alpha = 0.6, size = 3.5,
                add = "reg.line", conf.int = TRUE, 
                cor.coef = T, cor.method = "pearson",
                xlab = "Disinfectants Peak Area", ylab = "Read Counts") +
  theme(axis.text.y = element_text(size=21),
        axis.text.x = element_text(size=21),
        axis.title.y = element_text(size=21),
        axis.title.x = element_text(size=21),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(1, "mm"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(data_log_filtered, "log_QAC_bt", fill = "gray")
yplot <- ggdensity(data_log_filtered, "counts", fill = "gray")+rotate()
# Cleaning the plots
sp <- sp + rremove("legend")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
p1 <- insert_xaxis_grob(sp, xplot, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right") 
ggdraw(p2)

####################################################
# REFERENCE-BASED METABOLOMICS FOR SOURCE TRACKING # 
####################################################

# Code to run all the MS/MS spectra with fastMASST: https://github.com/robinschmid/microbe_masst/tree/v1.2.0
# Output for each feature is several tsv files, showing in which microbial species and food types the spectrum shows up
# Example tsv files are in the folder /MASST_pos_raw_data_example. The entire folder is too large to be uploaded to Github.
# The code below classified each feature into different source categories based on raw tsv files from the folder

path <- '/Users/ninazhao/Desktop/3DMM_MASST_raw/pos'
files <- list.files(path)
foods <- c()
microbes <- c()
homo <- c()
blankqc <- c()
pcp <- c()

for (file in files) {
  sp <- unlist(str_split(file, "_"))
  k <- sp[5]
  if (sp[6] == "counts") {
    df <- read.delim(file.path(path, file))
    source <- unlist(str_split(sp[7], "\\."))[1]
    if (source == "food") {
      foods <- c(foods, k) # Food as a possible source if the feature matched to foodMASST
    } else if (source == 'microbe') {
      if ('Homo sapiens' %in% df$Taxaname_file) {
        homo <- c(homo, k) # Host as a possible source if the feature matched to human cell lines in microbeMASST
        if (length(unique(df$Taxaname_file)) > 1) {
          microbes <- c(microbes, k) # Microbes as a possible source if the feature not only matched to human cell lines in microbeMASST
        }
      } else {
        microbes <- c(microbes, k) # Microbes as a possible source if the feature matched to microbeMASST
      }
      if ("qc" %in% df$Taxa_NCBI | "blank" %in% df$Taxa_NCBI) {
        blankqc <- c(blankqc, k) # May be process blanks if the feature matched to blank/qc in microbeMASST
      }
    }
  }
  else if (sp[6] == "matches.tsv"){
    df_match <- read.delim(file.path(path, file))
    if (any(grepl("MSV000080464|MSV000080473|MSV000080556|MSV000080557|MSV000082171|MSV000082173|MSV000081581", 
                  df_match$USI))){
      df_match_filtered <- df_match[grepl("MSV000080464|MSV000080473|MSV000080556|MSV000080557|MSV000082171|MSV000082173|MSV000081581", 
                                          df_match$USI),]
      if(!any(grepl("Blank", df_match_filtered$USI))){
        pcp <- c(pcp, k) 
        # Personal care products as possible sources if the feature is in PCP reference datasets and not in blank samples of these datasets
      }
    } 
  }
}

# Final classification of features
food_unique <- setdiff(foods, c(microbes, homo, blankqc, pcp)) # "Food" if the feature only matches to FoodMASST
microbe_unique <- setdiff(microbes, c(foods, homo, blankqc, pcp)) # "Microbe" if the feature only matches to MicrobeMASST
homo_unique <- setdiff(homo, c(foods, microbes, blankqc, pcp)) # "Human cell lines" if the feature only matches to human cell line samples in MicrobeMASST
pcp_unique <- pcp # "PCP" if the feature matches to PCP references
blankqc <- setdiff(blankqc, pcp) # "Blank/QC" if the feature matches to blank or qc samples in MicrobeMASST
multiple <- setdiff(unique(c(foods, microbes, homo)),
                    c(food_unique, microbe_unique, homo_unique, pcp_unique, blankqc)) # "Multiple" if the feature is not in blanks and matched to multiple categories

# Export tsv file for node labeling in Cytoscape
feature_each_class <- 
  data.frame(feature.ID = c(food_unique, microbe_unique, homo_unique, pcp_unique, blankqc, multiple),
             class = c(rep("Food", length(food_unique)), rep("Microbial", length(microbe_unique)),
                       rep("Human_cell_line", length(homo_unique)), rep("PCP", length(pcp_unique)),
                       rep("Blank_QC", length(blankqc)), rep("Multiple", length(multiple))))
write_tsv(feature_each_class, "MASST_analysis_results/3DMM_pos_reference_based_source_tracking.tsv")

# Source-tracking for module-specific features
module_count <- read.csv("Feature_tables/3DMM_pos_feature_module_count.csv")
module_count$feature.ID <- as.character(module_count$feature.ID)
module_count_filtered <- module_count %>% dplyr::filter(rowSums(.[, -1] > 0) == 1) %>%
  mutate(Module = colnames(.[, -1])[apply(.[, -1] > 0, 1, which)]) %>%
  left_join(., feature_each_class, by = "feature.ID") # Filter for module-specific feature and merge with MASST-based sources
module_count_filtered$class[is.na(module_count_filtered$class)] <- "Unknown"
summary_table <- as.data.frame(table(module_count_filtered$Module, module_count_filtered$class))
write.csv(summary_table, "MASST_analysis_results/3DMM_pos_module_specific_feature_source_tracking.csv")

########################################################
# FEATURE TABLE BLANK SUBTRACTION - NEGATIVE IONIZATION#
########################################################

metadata <- read.csv("Metadata/3DMM_metadata.csv")
feature_neg <- read.csv(file = "Feature_tables/3DMM_neg_quant.csv")
colnames(feature_neg) <- gsub(".mzML.Peak.area", "", colnames(feature_neg))
metadata_filtered <- metadata[metadata$Well %in% colnames(feature_neg),]

# For each plate, check if sample peak area is five times higher than the average of blanks and impute as zero if not.
plate <- c("P1_","P2_","P3_","P4_","P5_","P6_","P7_","P8_","P9_","P10_")
feature_neg_subtract <- feature_neg[,c("row.ID", "row.m.z", "row.retention.time")]
for (i in 1:length(plate)){
  metadata_plate <- metadata_filtered[grepl(plate[i], metadata_filtered$Well),]  # map for specific plate
  swabblk_plate <- metadata_plate$Well[metadata_plate$Barcode == "swab.blank"]  # name vector for swab blk in this plate
  sample_plate <- metadata_plate$Well[!metadata_plate$Barcode %in% c("swab.blank","solvent.blank")]  # name vector for sample in this plate
  swabblk_avr <- rowMeans(feature_neg[,swabblk_plate]) # average of swab blk
  feature_sample_plate <- feature_neg[,sample_plate] 
  t <- (feature_sample_plate > 5*swabblk_avr)  # if higher than 5*blks
  feature_sample_plate_subtracted <- t * feature_sample_plate - swabblk_avr  #subtract blk
  feature_sample_plate_subtracted <- (feature_sample_plate_subtracted > 0) * feature_sample_plate_subtracted
  feature_neg_subtract <- cbind(feature_neg_subtract, feature_sample_plate_subtracted)
}

# remove feature with less than two detections in all modules
feature_neg_filtered <- feature_neg_subtract %>% column_to_rownames("row.ID") %>% dplyr::select(contains("_")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID")
feature_metadata <- merge(feature_neg_filtered, metadata_filtered[,c("Well", "Module")], 
                          by.x = "SampleID", by.y = "Well", all.x = T, all.y = F)
feature_module_count <- feature_metadata %>% column_to_rownames("SampleID") %>% 
  aggregate(list(feature_metadata$Module), function(c)sum(c!=0)) %>% 
  dplyr::select(!contains("Module")) %>% column_to_rownames("Group.1")
feature_keep <- colnames(feature_module_count)[apply(feature_module_count, 2, function(x) any(x > 1))]

# save feature module count table for future processing
feature_module_count_filtered <- feature_module_count %>% 
  dplyr::select(all_of(feature_keep)) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("feature.ID")
write.csv(feature_module_count_filtered, "Feature_tables/3DMM_neg_feature_module_count.csv", row.names = F)

# save blank removed feature table for FBMN job on GNPS
feature_neg_duplicate_filtered <-  feature_neg_subtract %>% dplyr::filter(feature_neg_subtract$row.ID %in% feature_keep)
# modify the feature table format to be compatible with GNPS FBMN workflow
colnames(feature_neg_duplicate_filtered)[1:3] <- c("row ID", "row m/z", "row retention time")
colnames(feature_neg_duplicate_filtered)[4:ncol(feature_neg_duplicate_filtered)] <- 
  paste(colnames(feature_neg_duplicate_filtered)[4:ncol(feature_neg_duplicate_filtered)], ".mzML Peak area", sep = "")
# export feature table for FBMN
write.csv(feature_neg_duplicate_filtered, "Feature_tables/3DMM_neg_quant_blk_duplicate_filtered.csv", row.names = F)

#################################################################
# FILTER MGF FILE TO REMOVE BLANK FEATURES - NEGATIVE IONIZATION#
#################################################################

# "clean_mgf_featureID" was defined in early parts of this notebook to clean up mgf based on feature numbers
# modify mgf file to remove blank features
mgf_raw <- readLines("MGF/3DMM_neg_FBMN.mgf") 
mgf_modified <- clean_mgf_featureID(mgf_raw, feature_keep)
writeLines(mgf_modified, "MGF/3DMM_neg_FBMN_blk_rm.mgf") # mgf file for FBMN job on GNPS and fastMASST analysis

##########################################
# HEATMAPS FOR FEATURES UNIQUE TO NODE 3 #
##########################################

# Positive ionization features
pos_feature <- read.csv("Feature_tables/3DMM_pos_quant_blk_duplicate_filtered.csv")
metadata <- read.csv("Metadata/3DMM_metadata.csv")
pos_module_count <- read.csv("Feature_tables/3DMM_pos_feature_module_count.csv")
Node_3_specific <- pos_module_count %>% dplyr::filter(rowSums(.[, -1] > 0) == 1) %>%
  mutate(Module = colnames(.[, -1])[apply(.[, -1] > 0, 1, which)]) %>%
  dplyr::filter(.[,"Module"] == "Node.3") # Find features specific to Node 3
pos_feature_node3 <- pos_feature %>% dplyr::filter(pos_feature$row.ID %in% Node_3_specific$feature.ID) # Filter the feature table 
Node3_samples <- paste(metadata$Well[metadata$Module == "Node 3" & !(grepl("Control", metadata$Rack))], ".mzML.Peak.area", sep = "")
pos_feature_node3_cleanup <- pos_feature_node3 %>%
  column_to_rownames("row.ID") %>%
  dplyr::select(any_of(Node3_samples)) %>%
  t() %>% as.data.frame() 
# Normalize each cell by dividing by the maximum value in its column
feature_normalized <- as.data.frame(lapply(pos_feature_node3_cleanup, function(col) col / max(col)))
feature_normalized_cleanup <- feature_normalized %>%
  mutate(Sample = sub(".mzML.Peak.area", "", rownames(pos_feature_node3_cleanup))) %>%
  left_join(., metadata[, c("Location", "Well")], by = c("Sample" = "Well"))
# Select the maximum value at each location of Node 3
feature_normalized_max <- feature_normalized_cleanup %>% group_by(Location) %>% summarise(across(starts_with("X"), max))
data <- feature_normalized_max %>% column_to_rownames("Location") %>% as.matrix()
rownames(data) <- c("A4", "A5 (OGS)", "D2 (Cupola)","D4","D5","F2 (Hatch to PMM)",
                    "F4 (WHC)","F5 (T2)","O2/A2 (ARED)", "O4","O5" )
my_heatmap <- pheatmap(data, cluster_rows = T,
                       border_color = NA, 
                       clustering_method = "ward.D2",
                       fontsize = 16, 
                       angle_col = 45, 
                       show_colnames = F, show_rownames = T,
                       cellwidth = 3, cellheight = 20)

# Negative ionization features
neg_feature <- read.csv("Feature_tables/3DMM_neg_quant_blk_duplicate_filtered.csv")
metadata <- read.csv("Metadata/3DMM_metadata.csv")
neg_module_count <- read.csv("Feature_tables/3DMM_neg_feature_module_count.csv")
Node_3_specific <- neg_module_count %>% dplyr::filter(rowSums(.[, -1] > 0) == 1) %>%
  mutate(Module = colnames(.[, -1])[apply(.[, -1] > 0, 1, which)]) %>%
  dplyr::filter(.[,"Module"] == "Node.3") # Find features specific to Node 3
neg_feature_node3 <- neg_feature %>% dplyr::filter(neg_feature$row.ID %in% Node_3_specific$feature.ID) # Filter the feature table 
Node3_samples <- paste(metadata$Well[metadata$Module == "Node 3" & !(grepl("Control", metadata$Rack))], ".mzML.Peak.area", sep = "")
neg_feature_node3_cleanup <- neg_feature_node3 %>%
  column_to_rownames("row.ID") %>%
  dplyr::select(any_of(Node3_samples)) %>%
  t() %>% as.data.frame() 
# Normalize each cell by dividing by the maximum value in its column
feature_normalized <- as.data.frame(lapply(neg_feature_node3_cleanup, function(col) col / max(col)))
feature_normalized_cleanup <- feature_normalized %>%
  mutate(Sample = sub(".mzML.Peak.area", "", rownames(neg_feature_node3_cleanup))) %>%
  left_join(., metadata[, c("Location", "Well")], by = c("Sample" = "Well"))
# Select the maximum value at each location of Node 3
feature_normalized_max <-  feature_normalized_cleanup %>% group_by(Location) %>% summarise(across(starts_with("X"), max))
order <- c("NOD3O2 (ARED) / NOD3A2 (ARED BAR PLATFORM OVER BEAM HATCH PORTAL)",
           "NOD3A5", "NOD3F4", "NOD3A4", "NOD3O5", "NOD3D5",
           "NOD3O4", "NOD3F2 (Hatch to PMM)", "NOD3F5 (T2)",
           "NOD3D2 (Cupola / Hatch to Cupola)", "NOD3D4")
feature_normalized_max <- feature_normalized_max[order(factor(feature_normalized_max$Location, levels = order)), ]
data <- feature_normalized_max %>% column_to_rownames("Location") %>% as.matrix()
rownames(data) <- c("A4", "A5 (OGS)", "D2 (Cupola)","D4","D5","F2 (Hatch to PMM)",
                    "F4 (WHC)","F5 (T2)","O2/A2 (ARED)", "O4","O5" )
my_heatmap <- pheatmap(data, cluster_rows = F,
                       border_color = NA, clustering_method = "ward.D2",
                       fontsize = 16, angle_col = 45, show_colnames = F, show_rownames = T,
                       cellwidth = 8, cellheight = 20)

#####################################################
# CO-MOLECULAR NETWORKING WITH URBANIZATION DATASET #
#####################################################
cluster <- read.csv("Co_network_urbanization/3DMM_urbanization_coMN_cluster_info.csv")
cluster_subset <- cluster[(cluster$G1 > 1) & (cluster$G3 > 1) & (cluster$G2 == 0),] #Group1: 3DMM sample; Group2: 3DMM blank; Group3: Urbanization sample
family_to_keep <- cluster_subset$componentindex[cluster_subset$componentindex != -1] 
cluster_to_keep <- cluster_subset$cluster.index
keep_test <- (cluster$cluster.index %in% cluster_to_keep) | (cluster$componentindex %in% family_to_keep) & (cluster$G2 == 0)
cluster_results <- cbind(cluster, keep_test)
metadata <- read.csv("Co_network_urbanization/3DMM_urbanization_coMN_metadata.csv")
metadata$subclass[metadata$dataset == "urbanization" & metadata$type == "blank"] <- "Urbanization-Blank"
metadata$subclass[metadata$dataset == "3DMM" & metadata$type == "blank"] <- "3DMM-Blank"
abundance <- read.csv("Co_network_urbanization/3DMM_urbanization_coMN_cluster_abundance.csv")
abundance_filtered <- abundance %>% column_to_rownames("X.OTU.ID") %>% t() %>% as.data.frame() %>% rownames_to_column("Sample")
abundance_merge <- merge(metadata, abundance_filtered, by.x = "file", by.y = "Sample")
abundance_merge <- abundance_merge[,!colnames(abundance_merge) %in% c("file","dataset","type")]
result <- abundance_merge %>%
  group_by(subclass) %>%
  summarise_all(sum)
result_filtered <- result %>% column_to_rownames("subclass") %>% t() %>% as.data.frame() %>% rownames_to_column ("cluster.index")
cluster_final <- cluster_results[,c("cluster.index", "keep_test")]
cluster_final <- merge(cluster_final, result_filtered, by = "cluster.index")
write.csv(cluster_final,"Co_network_urbanization/Common_cluster_group_abundance.csv") # file for Cytoscape visualization

#############################################################
# CORRELATION NETWORK BETWEEN METABOLOMICS AND 16S FEATURES #
#############################################################
metadata <- read.csv("Metadata/3DMM_metadata.csv")
feature_pos <- read.csv(file = "Feature_tables/3DMM_pos_quant.csv")
colnames(feature_pos) <- gsub(".mzML.Peak.area", "", colnames(feature_pos))
metadata_filtered <- metadata[metadata$Well %in% colnames(feature_pos),]

# perform blank subtraction on metabolomics feature table without imputation to better fit correlation analysis
feature_filter <- feature_pos %>% column_to_rownames("row.ID") %>% dplyr::select(contains("_"))
feature_blank <- dplyr::select(feature_filter, any_of(metadata_filtered$Well[metadata_filtered$Rack == "Field Control"]))
feature_sample <- dplyr::select(feature_filter, any_of(metadata_filtered$Well[metadata_filtered$Rack != "Negative Control" & 
                                                                                metadata_filtered$Rack != "Field Control" &
                                                                                metadata_filtered$Rack != "solvent.blank" &
                                                                                metadata_filtered$Rack != "swab.blank"]))
all(rownames(feature_blank) == rownames(feature_sample))
blk_check <- data.frame(avg_blank=rowMeans(feature_blank), avg_sample=rowMeans(feature_sample), LOD=3*rowMeans(feature_blank))
blk_check$ratio_sample_blank <- (blk_check$avg_sample+1)/(blk_check$avg_blank+1)
blk_check$bin <- ifelse(blk_check$ratio_sample_blank > 3, 1, 0)
all(rownames(blk_check) == rownames(feature_sample))
feature_blk_removed <- feature_sample[blk_check$bin == 1,]
feature_filter <- feature_blk_removed %>% t() %>% as.data.frame() %>% rownames_to_column("Well") %>%
  left_join(metadata %>% dplyr::select(c("Well", "Labels")), by = "Well") %>% 
  mutate(Labels = gsub("[-/]", ".", Labels)) %>%
  column_to_rownames("Labels") %>% dplyr::select(-Well)
colnames(feature_filter) <- paste("F", colnames(feature_filter), sep = "")

# Non-zero variance filter and RCLR transformation
feature_filter_nzv <- feature_filter %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE))))
feature_rclr <- decostand(feature_filter_nzv, method = "rclr") %>% rownames_to_column("ID")
# Non-zero variance filter and RCLR transformation for 16S features
microbiome <- read.csv("16s_data/katharo_deblur_gg2_unrarefied_genera_collapsed.csv", row.names = 1)
rownames(microbiome) <- sub("151322.14542.", "", rownames(microbiome))
rownames(microbiome) <- sub("150297.14542.", "", rownames(microbiome))
microbiome_nzv <- microbiome %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE))))
microbiome_rclr <- decostand(microbiome_nzv, method = "rclr") %>% rownames_to_column("ID")

pkgs <- c("tidyverse", "ggpubr","santaR",  "mixOmics", "plotly", "RColorBrewer", "magrittr", "Hmisc", "Rmisc",
          "ALDEx2", "phyloseq", "CoDaSeq", "reshape2", "vegan", "microbiome", "DivNet", "igraph", "ggraph")
for (pkg in pkgs) {library(pkg, character.only = TRUE)}

# Order based on sample ID
metaboliteINT <- feature_rclr[order(feature_rclr$ID),]
microbiomeINT <- microbiome_rclr[order(microbiome_rclr$ID),]
# Filter
tokeep <- inner_join(metaboliteINT, microbiomeINT, by = "ID") %>% dplyr::select(ID)
metaboliteINT_filter <- dplyr::filter(metaboliteINT, ID %in% tokeep$ID)
microbiomeINT_filter <- dplyr::filter(microbiomeINT, ID %in% tokeep$ID)

# Check that filter data are in the same order
all(metaboliteINT_filter$ID == microbiomeINT_filter$ID)

# read in microbeMASST results
microbeMASST <- read.csv("MASST_analysis_results/3DMM_pos_reference_based_source_tracking.tsv", sep = "\t")
microbe_feature <- microbeMASST$feature.ID[microbeMASST$class == "Microbial"]
microbe_feature <- paste("F", microbe_feature, sep = "")
sum(microbeMASST$feature.ID[microbeMASST$class == "Microbial"] %in% rownames(feature_blk_removed))

# Remove ID column
metaboliteINT_final <- metaboliteINT_filter %>% dplyr::select(-c("ID")) %>% dplyr::select(any_of(microbe_feature))
microbiomeINT_final <- microbiomeINT_filter %>% select(-c("ID"))

# extract 16S taxanomic information
extract_value <- function(input_string, prefix) {
  pattern <- paste0(prefix, "([^\\.]+)")
  match <- str_match(input_string, pattern)
  return(match[2])
}
prefixes <- c("d__", "p__", "c__", "o__", "f__", "g__")
process_string <- function(input_string) {
  values <- sapply(prefixes, extract_value, input_string = input_string)
  return(as.list(values))
}
processed_data <- lapply(colnames(microbiomeINT_final), process_string)
microbiome_metadata <- bind_rows(processed_data)
colnames(microbiome_metadata) <- c("d", "p", "c", "o", "f", "g")
microbiome_metadata <- microbiome_metadata %>% mutate(featureID = paste0("M", 1:nrow(microbiome_metadata)))
colnames(microbiomeINT_final) <- microbiome_metadata$featureID
sort(unique(microbiome_metadata$p))
microbiome_metadata$p[grepl("Firmicutes", microbiome_metadata$p)] <- "Firmicutes"
microbiome_metadata$p[grepl("Myxococcota", microbiome_metadata$p)] <- "Myxococcota"

# Generate matrices for pair-wise correlation between the datasets
matrix_MetMic <- cbind(metaboliteINT_final, microbiomeINT_final)
# Calculate correlations
corr_MetMic <- rcorr(as.matrix(matrix_MetMic), type = "pearson")

# Function - flatten correlation matrix (from sthda website)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
# Flatten correlation matrices, adjust with BH and filter for FDR < 0.05
corr_MetMic_flat <- flattenCorrMatrix(corr_MetMic$r,corr_MetMic$P)
corr_MetMic_flat$padj <- p.adjust(corr_MetMic_flat$p, method = "BH")
corr_MetMic_final <- dplyr::filter(corr_MetMic_flat, 
                                   corr_MetMic_flat$padj < 0.05)

# Filter the network to keep only positive correlations
corr_MetMic_final_filter <- corr_MetMic_final %>%
  filter(grepl("F", corr_MetMic_final$row) &
           grepl("M", corr_MetMic_final$column) &
           cor >0)
result <- corr_MetMic_final_filter %>% left_join(microbiome_metadata, by = c("column" = "featureID"))
# Generate nodes
met_nodes <- colnames(metaboliteINT_final)
bac_nodes <- colnames(microbiomeINT_final)
# Assign category at each node
met_nodes <- data.frame(node = met_nodes, category = rep("Metabolite",length(met_nodes)), phylum = rep("Metabolite",length(met_nodes)))
bac_nodes <- data.frame(node = bac_nodes, category = rep("Bacteria",length(bac_nodes)))
bac_nodes <- bac_nodes %>% left_join(microbiome_metadata %>% dplyr::select(c("p", "featureID")), by = c("node" = "featureID"))
colnames(bac_nodes)[3] <- "phylum"

met_nodes_filter <- met_nodes %>% dplyr::filter(node %in% corr_MetMic_final_filter$row)
bac_nodes_filter <- bac_nodes %>% dplyr::filter(node %in% corr_MetMic_final_filter$column)
nodes <- rbind(met_nodes_filter, bac_nodes_filter)
node_colors <- c("#83A4C9", "#5F8AB9")
names(node_colors) <- c("Metabolite", "Bacteria")
# Generate network table
network_table <- corr_MetMic_final_filter
network_table <- network_table[complete.cases(network_table),]
colnames(network_table)[1:3] <- c("from", "to", "correlation")
network_table$direction <- ifelse(network_table$correlation < 0, "neg", "pos")
network_table$weight <- abs(network_table$correlation)
# Generate network
network_whole <- graph_from_data_frame(d = network_table, vertices = nodes, directed = FALSE)
phylum_colors <- brewer.pal(n = length(unique(V(network_whole_pos)$phylum)), name = "Paired")

network_pos_3DMM <- ggraph(network_whole, layout = "fr") +
  geom_edge_link(aes(color = correlation), edge_width = 0.5) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-0.35, 0.35)) +
  geom_node_point(aes(shape = category, col = phylum), size = 3) +
  geom_node_text(aes(label = name), size = 3, repel = TRUE) +
  # scale_color_manual(values = phylum_colors) +
  theme_graph(base_family = "Helvetica", title_size = 10) +
  labs(title = "Microbiome-Disinfectants Network",
       shape = "Category", edge_color = "Pearson Correlation")

# Modify network on Cytoscape
library(RCy3)
cytoscapePing ()
cytoscapeVersionInfo ()
createNetworkFromIgraph(network_whole, "network_pos_3DMM")
