##############################################################################################
# Script Name:   OPLS_VIP_Analysis.R
# Description:   Performs OPLS-DA based VIP ranking of proteins for
#                multi-class proteomics data on sCJD subtypes.  
#                Performs CJD vs CTRL, subtype to subtype and subtype vs rest comparisons. 
#                Generates VIP score barplots for the top proteins,
#                highlighting proteins previously selected by RF analyses.
#
# Requirements:  R packages: tidyverse, readxl, ggtext, ropls, rstudioapi
#
# Input:         Dataframe 'df' with protein abundances (columns)
#                and grouping variables 'Group' and 'SubGroup' (factor or character).
#
# Output:        Barplot PNG files saved in a 'results' directory next to this script.
#                CSV files containing model scores and VIP score rankings, stored in 'results'
#
# Usage:         Run this script in RStudio or source it after loading data.
#
# Author:        Isabel Houtkamp
# Date:          2025-08-08
###############################################################################################


# Load required packages

library(tidyverse)
library(readxl)
library(ropls)
library(ggtext)
library(rstudioapi)

# set seed
set.seed(114)

# Extract script directory a
script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

# Define results folder path
results_dir <- file.path(script_dir, "opls_results")

# Create 'results' folder if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Load curated data
df <- read_excel(paste0(script_dir, "/sCJD-subtypes-main/data/curated/olink.xlsx"))

# Specify columns to drop from full dataframe (sample info)
columns_to_drop = c('SampleID', 'age at LP', 'Sex', 'Codon 129',
                    'onset-LP', 'onset-death', 'LP-death', 'Group', 'Strain', 'NP_subtype', 'SubGroup', 
                    'group_binary')


# List protein names of top20 protein panel resulting from RF analyses
highlight_proteins <- c("FOSB", "PSIP1", "HEXIM1", "PARP-1", "APEX1", "MAPT", "NEFL", "WASF1", "ARHGEF12", "HDGF", "PAG1", "GPC5", "CAMKK1", "PPP3R1", "FKBP4", "EIF4G1", "THOP1", "METAP1D", "PRDX3", "CCDC80")

### CJD (all subtypes) vs controls comparison ###

# Create input protein data (X) and grouping vector (Y)
X = df %>% select(-any_of(columns_to_drop))
Y = as.factor(df$Group)

# Fit OPLS-DA (binary)
opls_model <- opls(X, Y, predI = NA, orthoI = NA, crossvalI = length(Y))

# Extract and rank VIP scores
vip_scores <- getVipVn(opls_model)
vip_ranking <- sort(vip_scores, decreasing = TRUE)

# Top 20 for plotting
top20 <- head(vip_ranking, 20)
df_plot <- data.frame(
  protein = names(top20),
  vip = as.numeric(top20),
  stringsAsFactors = FALSE
)


# Order for horizontal bars
df_plot <- df_plot[order(df_plot$vip), ]
df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)

# Annotate panel proteins
df_plot$highlight <- df_plot$protein %in% highlight_proteins

# Plot
p <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
  theme_classic() + 
  theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
  labs(
    x = "VIP score",
    title = paste("Top 20" , "VIP Scores for sCJD vs CTRL")
  )


# Save plot
ggsave(filename = file.path(results_dir, "Top20_VIP_CJD_vs_CTRL.png"), 
       plot = p, width = 7, height = 5, dpi = 1200)

# Repeat with top 10 only

# Top 20 for plotting
top10 <- head(vip_ranking, 10)
df_plot <- data.frame(
  protein = names(top10),
  vip = as.numeric(top10),
  stringsAsFactors = FALSE
)

# Order for horizontal bars
df_plot <- df_plot[order(df_plot$vip), ]
df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)

# Annotate panel proteins
df_plot$highlight <- df_plot$protein %in% highlight_proteins


p10 <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
  theme_classic() + 
  theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
  labs(
    x = "VIP score",
    title = paste("Top 10", "VIP Scores for sCJD vs CTRL")
  )

# Save plot
ggsave(filename = file.path(results_dir, "Top10_VIP_CJD_vs_CTRL.png"), 
       plot = p10, width = 7, height = 5, dpi = 1200)


## Save model quality
quality_df <- as.data.frame(t(opls_model@summaryDF))
write.csv(quality_df,
          file = file.path(results_dir, paste0("ModelQuality_CJDvsCTRL.csv")),
          row.names = TRUE)

## Save VIP ranking 
vip_scores <- getVipVn(opls_model)
vip_ranking <- sort(vip_scores, decreasing = TRUE)
vip_df <- data.frame(Protein = names(vip_ranking), VIP = vip_ranking)
write.csv(vip_df,
          file = file.path(results_dir, paste0("VIP_Ranking_CJDvsCTRL.csv")),
          row.names = FALSE)


#### All to all subgroup comparisons ###

# list subgroups and all possible pairings
subgroups <- unique(df$SubGroup)
pairs <- combn(subgroups, 2, simplify = FALSE)  # all subgroup pairs

for (pair in pairs) {
  
  group1 <- pair[1]
  group2 <- pair[2]
  
  # Filter only the two subgroups for comparison
  df_pair <- df %>% filter(SubGroup %in% c(group1, group2))
  
  X = df_pair %>% select(-any_of(columns_to_drop))
  Y = as.factor(df_pair$SubGroup)
  
  # Fit OPLS-DA (binary)
  opls_model <- opls(X, Y, predI = NA, orthoI = NA, crossvalI = length(Y))
  
  # Extract and rank VIP scores
  vip_scores <- getVipVn(opls_model)
  vip_ranking <- sort(vip_scores, decreasing = TRUE)
  
  # Top 20 for plotting
  top20 <- head(vip_ranking, 20)
  df_plot <- data.frame(
    protein = names(top20),
    vip = as.numeric(top20),
    stringsAsFactors = FALSE
  )
  
  # Order for horizontal bars
  df_plot <- df_plot[order(df_plot$vip), ]
  df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)
  
  # Annotate panel proteins
  df_plot$highlight <- df_plot$protein %in% highlight_proteins
  
  # Plot
  p <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
    theme_classic() + 
    theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
    labs(
      x = "VIP score",
      title = paste("Top 20" , "VIP Scores for", as.character(group1), "vs", as.character(group2))
    )
  
  # Top 10 for plotting
  top10 <- head(vip_ranking, 10)
  df_plot <- data.frame(
    protein = names(top10),
    vip = as.numeric(top10),
    stringsAsFactors = FALSE
  )
  
  # Order for horizontal bars
  df_plot <- df_plot[order(df_plot$vip), ]
  df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)
  
  # Annotate panel proteins
  df_plot$highlight <- df_plot$protein %in% highlight_proteins
  
  # Plot
  p10 <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
    theme_classic() + 
    theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
    labs(
      x = "VIP score",
      title = paste("Top 10" , "VIP Scores for", as.character(group1), "vs", as.character(group2))
    )
  
  print(p)  # show plot
  print(p10)
  
  # save each plot
  ggsave(filename = file.path(results_dir, paste0("VIP_", group1, "_vs_", group2, ".png")), plot = p, width = 7, height = 5, dpi = 1200)
  ggsave(filename = file.path(results_dir, paste0("VIP_", group1, "_vs_", group2, "_top10.png")), plot = p10, width = 7, height = 5, dpi = 1200)
  
  
  ## Save model quality
  quality_df <- as.data.frame(t(opls_model@summaryDF))
  write.csv(quality_df,
            file = file.path(results_dir, paste0("ModelQuality_", group1, "_vs_", group2, ".csv")),
            row.names = TRUE)
  
  ## Save VIP ranking 
  vip_scores <- getVipVn(opls_model)
  vip_ranking <- sort(vip_scores, decreasing = TRUE)
  vip_df <- data.frame(Protein = names(vip_ranking), VIP = vip_ranking)
  write.csv(vip_df,
            file = file.path(results_dir, paste0("VIPranking", group1, "_vs_", group2, ".csv")),
            row.names = FALSE)
  
}


### Subgroup vs Rest comparisons ###

# Keep subtype data only (drop CTRLS)
df_sub <- df %>% filter(SubGroup != "CTRL")
subgroups <- unique(df_sub$SubGroup)

# Loop over all subgroups for a subgroup vs Rest comparison. 
for (sub in subgroups) {
  
  # Make binary group: current subgroup vs "Rest"
  df_pair <- df_sub
  df_pair$group_binary <- ifelse(df_pair$SubGroup == sub, sub, "Rest")
  
  X = df_pair %>% select(-any_of(columns_to_drop))
  Y = as.factor(df_pair$group_binary)
  
  # Fit OPLS-DA
  opls_model <- opls(X, Y, predI = NA, orthoI = NA, crossvalI = length(Y))
  
  # VIP ranking
  vip_scores <- getVipVn(opls_model)
  vip_ranking <- sort(vip_scores, decreasing = TRUE)
  
  # Top 20 for plotting
  top20 <- head(vip_ranking, 20)
  df_plot <- data.frame(
    protein = names(top20),
    vip = as.numeric(top20),
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_plot[order(df_plot$vip), ]
  df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)
  df_plot$highlight <- df_plot$protein %in% highlight_proteins
  
  # Plot
  p <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
    labs(
      x = "VIP score",
      title = paste("Top 20 VIP Scores for", sub, "vs Rest")
    )
  
  # Repeat for top 10
  top10 <- head(vip_ranking, 10)
  df_plot <- data.frame(
    protein = names(top10),
    vip = as.numeric(top10),
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_plot[order(df_plot$vip), ]
  df_plot$protein <- factor(df_plot$protein, levels = df_plot$protein)
  df_plot$highlight <- df_plot$protein %in% highlight_proteins
  
  # Plot
  p10 <- ggplot(df_plot, aes(x = vip, y = protein, fill = highlight)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "darkorange", "FALSE" = "dodgerblue4")) +
    theme_classic() +
    theme(axis.title.y = element_blank()) + xlim(c(0,5)) +
    labs(
      x = "VIP score",
      title = paste("Top 10 VIP Scores for", sub, "vs Rest")
    )
  
  # Save
  ggsave(filename = file.path(results_dir, paste0("VIP_", sub, "_vs_Rest.png")),
         plot = p, width = 7, height = 5, dpi = 1200)
  
  ggsave(filename = file.path(results_dir, paste0("VIP_", sub, "_vs_Rest_top10.png")),
         plot = p10, width = 7, height = 5, dpi = 1200)
  
  
  ## Save model quality
  quality_df <- as.data.frame(t(opls_model@summaryDF))
  write.csv(quality_df,
            file = file.path(results_dir, paste0("ModelQuality_", sub, "_vs_Rest.csv")),
            row.names = TRUE)
  
  ## Save VIP ranking 
  vip_scores <- getVipVn(opls_model)
  vip_ranking <- sort(vip_scores, decreasing = TRUE)
  vip_df <- data.frame(Protein = names(vip_ranking), VIP = vip_ranking)
  write.csv(vip_df,
            file = file.path(results_dir,paste0("VIP_Ranking", sub, "_vs_Rest.csv")),
            row.names = FALSE)
  
}

################################################################################
# End of Script
# Notes:
#   - Check 'opls_results' folder for generated PNG plots, model quality metrics, 
#     and VIP rankings of all proteins. 
#   - Ensure that 'df' is loaded and formatted correctly before running.
################################################################################

