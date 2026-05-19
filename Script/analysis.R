# ==============================================================================
#      RNA-Seq analysis
#      Infection of THP-1 macrophages with Helicobacter pylori B8 strain
# ==============================================================================
# Author: Dr. Angelika Lahnsteiner
# Date: 2025-11-22
# 
# ==============================================================================


# ------------------------------
# Load libraries
# ------------------------------
library("DESeq2")
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(eulerr)
library(grid)
library(dplyr)
library(tidyr)
library(pheatmap)
library(DOSE)
library(pathview)
library(VennDiagram)
library(ggvenn)
library(enrichplot)
library(circlize)


# ------------------------------
# prepare tables
# ------------------------------

# first load all count files from feature counts output from Galaxy
files <- list.files(pattern = "\\.tabular$", full.names = TRUE)

for (f in files) {
  name <- tools::file_path_sans_ext(basename(f))
  
  # Remove prefix
  name <- sub("^featureCounts on ", "", name)
  
  # Keep only up to 6h or 24h
  name <- sub("((6h|24h)).*", "\\1", name)
  
  assign(name, read.delim(f), envir = .GlobalEnv)
}


# rename colnames of each file for all loaded files
# Get all data frame names
df_names <- ls()

for (df in df_names) {
  
  # Check if object is a data frame
  if (is.data.frame(get(df))) {
    
    # Get the data frame
    temp <- get(df)
    
    # Rename first two columns
    colnames(temp)[1:2] <- c("GeneID", df)
    
    # Save back to environment
    assign(df, temp, envir = .GlobalEnv)
  }
}


#merge count files
df_list <- list(WT_EI_6h, WT_EII_6h, WT_EIII_6h, 
                WT_EI_24h, WT_EII_24h, WT_EIII_24h,
               Mock_EI_6h, Mock_EII_6h, Mock_EIII_6h,
               Mock_EI_24h, Mock_EII_24h, Mock_EIII_24h) #


counts <- Reduce(function(x, y) merge(x, y, by = "GeneID"), df_list)

# remove unused files to clean up environment
rm(df,df_names,f,files,name,temp, df_list,WT_EI_24h, WT_EI_6h, WT_EII_24h, WT_EII_6h, WT_EIII_24h, WT_EIII_6h,
   Mock_EI_24h, Mock_EI_6h, Mock_EII_24h, Mock_EII_6h, Mock_EIII_24h, Mock_EIII_6h)

  
# save the counts table for later usage
write.table(counts, file = "merged_counts_files.tsv", row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)

# remove col. GeneID and use it as rownames
rownames(counts) <- counts$GeneID
counts$GeneID <- NULL

# all columns ending in "_6h"
cols_6h <- grep("_6h$", colnames(counts), value = TRUE)
counts_6h <- counts[ , cols_6h ]

# all columns ending in "_24h"
cols_24h <- grep("_24h$", colnames(counts), value = TRUE)
counts_24h <- counts[ , cols_24h ]

# prepare a sample table
samples <- colnames(counts)

coldata <- data.frame(
  sample = samples,
  
  condition = ifelse(grepl("^Mock", samples), "Mock", "WT"),
  
  replicate = sub(".*_(EI+)_.*", "\\1", samples),
  
  time = sub(".*_(6h|24h).*", "\\1", samples)
)

# final sample table to build the contrasts
coldata

# ==============================================================================
# 1. Analysis of conditions across time points
# ==============================================================================

### 1) Split counts by timepoint using colnames --------------------------------
# columns ending with "_6h"
cols_6h   <- grep("_6h$",   colnames(counts), value = TRUE)
# columns ending with "_24h"
cols_24h  <- grep("_24h$",  colnames(counts), value = TRUE)
# extract the counts for _6h and _24h samples
counts_6h  <- counts[ , cols_6h ]
counts_24h <- counts[ , cols_24h ]

# subset coldata accordingly
cols_6h  <- coldata$time == "6h"
cols_24h <- coldata$time == "24h"
# extract the sample information for _6h and _24h samples
coldata_6h  <- droplevels(coldata[cols_6h, ])
coldata_24h <- droplevels(coldata[cols_24h, ])

### 2) Set Mock as reference for each timepoint --------------------------------
coldata_6h$condition  <- relevel(factor(coldata_6h$condition),  ref = "Mock")
coldata_24h$condition <- relevel(factor(coldata_24h$condition), ref = "Mock")

### 3) Build DESeq2 objects with design = ~ condition -----------------
dds_6h <- DESeqDataSetFromMatrix(
  countData = counts_6h,
  colData   = coldata_6h,
  design    = ~ condition
)

dds_24h <- DESeqDataSetFromMatrix(
  countData = counts_24h,
  colData   = coldata_24h,
  design    = ~ condition
)

### 4) Run DESeq2 for each timepoint -------------------------------------------
dds_6h  <- DESeq(dds_6h)
dds_24h <- DESeq(dds_24h)

### 5) Extract pairwise contrasts within each timepoint ------------------------
# ---- 6h ----
WT_vs_Mock_6h   <- results(dds_6h,  contrast = c("condition", "WT",  "Mock"))
# ---- 24h ----
WT_vs_Mock_24h  <- results(dds_24h, contrast = c("condition", "WT",  "Mock"))

### 6) Helper function to annotate Ensembl IDs, Gene symbols and save as TSV ---
annotate_and_save <- function(res, prefix) {
  res_df <- as.data.frame(res)
  
  # clean Ensembl IDs (strip version) from rownames
  ensembl_clean <- sub("\\..*", "", rownames(res_df))
  
  symbols <- mapIds(
    org.Hs.eg.db,
    keys      = ensembl_clean,
    column    = "SYMBOL",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  
  entrez <- mapIds(
    org.Hs.eg.db,
    keys      = ensembl_clean,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  
  res_annot <- data.frame(
    Ensembl = ensembl_clean,
    Symbol  = symbols,
    Entrez  = entrez,
    res_df,
    row.names = NULL
  )
  
  # sort by adjusted p-value if present
  if ("padj" %in% colnames(res_annot)) {
    res_annot <- res_annot[order(res_annot$padj), ]
  }
  
  ## a) Save full annotated table
  fname_full <- paste0(prefix, ".tsv")
  write.table(
    res_annot,
    file      = fname_full,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  ## b) Save filtered tables
  if (all(c("log2FoldChange", "padj") %in% colnames(res_annot))) {
    
    # upregulated
    res_up <- subset(res_annot,
                     log2FoldChange > 1 & !is.na(padj) & padj < 0.05)
    
    # downregulated
    res_down <- subset(res_annot,
                       log2FoldChange < -1 & !is.na(padj) & padj < 0.05)
    
    # combined significant
    res_sign <- subset(res_annot,
                       abs(log2FoldChange) > 1 & !is.na(padj) & padj < 0.05)
    
    fname_sign <- paste0(prefix, ".sign.tsv")
    fname_up   <- paste0(prefix, ".sign.up.tsv")
    fname_down <- paste0(prefix, ".sign.down.tsv")
    
    write.table(res_sign, fname_sign,
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(res_up, fname_up,
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(res_down, fname_down,
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  message("Saved: ",
          fname_full, ", ",
          if (exists("fname_sign")) fname_sign else "",
          if (exists("fname_up"))   paste0(", ", fname_up)   else "",
          if (exists("fname_down")) paste0(", ", fname_down) else "")
  
  invisible(res_annot)
}


##  7) Annotate & save all contrasts -------------------------------------------
# ---- 6h ----
annotate_and_save(WT_vs_Mock_6h, "DESeq2_WT_vs_Mock_6h_annotated")

# ---- 24h ----
annotate_and_save(WT_vs_Mock_24h, "DESeq2_WT_vs_Mock_24h_annotated")


##  8) Identify shared and uniquely differentially regulated genes -------------

# load DEGs
sign6h <- read.delim("~/Desktop/Analysis H.pylori/Datasets/DESeq2_WT_vs_Mock_6h_annotated.sign.tsv")
sign24h <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.sign.tsv")

# find shared/unique transcripts
shared <-  sign6h %>%
  filter(Ensembl %in% sign24h$Ensembl)

unique6h <- sign6h %>%
  filter(!Ensembl %in% sign24h$Ensembl)

unique24h <-  sign24h %>%
  filter(!Ensembl %in% sign6h$Ensembl)

# write tables
write.table(shared, "~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_annotated_sign_shared_6h_24h.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique6h, "~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_annotated_sign_unique_6h.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(unique24h, "~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_annotated_sign_unique_24h.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Plot them as Venn Diagram 

# --- Build exact sets so areas are correct ---
set_6h  <- unique(c(shared$Ensembl, unique6h$Ensembl))
set_24h <- unique(c(shared$Ensembl, unique24h$Ensembl))

# Quick sanity counts
n_unique6  <- sum(!set_6h %in% set_24h)    # exclusive to 6h
n_unique24 <- sum(!set_24h %in% set_6h)    # exclusive to 24h
n_shared   <- length(intersect(set_6h, set_24h))
cat("counts: unique6 =", n_unique6, " unique24 =", n_unique24, " shared =", n_shared, "\n")


# --- Fit area-proportional two-set diagram ---
fit <- euler(list("6h" = set_6h, "24h" = set_24h))

# --- Aesthetic parameters for a manuscript-quality figure ---
fills <- list(fill = c("#2C7FB8", "#A6CEE3"), alpha = 0.6)  # blue / orange, semi-transparent
labels <- list(
  font = 2,              # bold labels
  labels = c("unique for 6h", "unique for 24h"),
  vjust = 0.5,
  cex = 1.1
)
quant_labels <- TRUE     # show numbers inside the diagram

# Plot the Venn Diagram
plot(fit,
     fills = fills,
     labels = labels,
     quantities = list(type = "counts", fontfamily = "sans", fontsize = 12),
     #main = "Differential genes: 6h vs 24h",
     main.cex = 1.1
)



# ==============================================================================
# 2. Heatmap
# ==============================================================================


#  Make sure colddata rows are ordered to match count columns
md <- coldata[match(colnames(counts), coldata$sample), ]

# 3) Create grouping variable (condition_time)
group_time <- with(md, paste(condition, time, sep = "_"))

# 4) Compute mean across replicates for each Group_Time
counts_mean <- sapply(
  split(seq_along(colnames(counts)), group_time),
  function(idx) rowMeans(counts[, idx, drop = FALSE])
)

rownames(counts_mean) <- sub("\\..*$", "", rownames(counts_mean))

# load transcripts based in 24h WT versus Mock 
DESeq2_WT_vs_Mock_24h_annotated.sign <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.sign.tsv")
sign_genes <-DESeq2_WT_vs_Mock_24h_annotated.sign[1:53,]

#remove a row with NAs in the symbol column
sign_genes <- na.omit(sign_genes) 
genes <- sign_genes$Ensembl
filtered_targets <- counts_mean[rownames(counts_mean) %in% genes,]
filtered_targets <- as.data.frame(filtered_targets)

#gene annotation
annotation <- sign_genes[c(1,2)]
colnames(annotation)[1] <- "GeneID" 
filtered_targets$GeneID <- rownames(filtered_targets)
filtered_targets <- merge(annotation,filtered_targets,by="GeneID")
genes <- filtered_targets$Symbol
filtered_targets <- filtered_targets[-c(1:2)]
rownames(filtered_targets) <- genes

# plot order (since we are not clustering columns)
plot_order <- c("Mock_6h", "Mock_24h", "WT_6h", "WT_24h")
filtered_targets <- filtered_targets[, plot_order]

# Annotation
annotation_row <- sign_genes[c(1,2)]
geneid <- annotation_row$Ensembl
annotation_row <- annotation_row[-c(1)]
rownames(annotation_row) <- geneid

# extract Time (6h / 24h) from column names
annotation_col <- data.frame(
  Time = sub(".*_(6h|24h).*", "\\1", colnames(counts_mean)),
  row.names = colnames(counts_mean),
  stringsAsFactors = FALSE
)
annotation_col <- annotation_col[plot_order, , drop = FALSE]

# change it to factor
annotation_col$Time <- factor(annotation_col$Time, levels = c("6h", "24h"))

# define colors
ann_colors <- list(
  Time = c("6h" = "#56B4E9", "24h" = "#E69F00")
)

color <- colorRampPalette((c( "darkblue","white", "red")))(50)

# plot heatmap
pheatmap(filtered_targets,
         cluster_rows = TRUE,         # Cluster rows (CpG sites)
         cluster_cols = FALSE,        # Keep columns in specified order
         scale = "row",               # Scale by row
         #annotation_row = annotation_row,  # Column annotation
         #main = "Top 50 DEGs in WT versus Mock based on 24h infection",
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = NA,
         color = color,
         annotation_col = annotation_col, 
         annotation_colors = ann_colors)



# ==============================================================================
# 3. GO enrichment analysis
# ==============================================================================

# Perform the analysis for both datasets 6h and 24h independently 
#----------
## Diff.regulated genes after 6h infection
genes.up <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_6h_annotated.sign.up.tsv")
genes.down <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_6h_annotated.sign.down.tsv")
allOE_genes <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_6h_annotated.tsv")

## Diff.regulated genes after 24h infection
genes.up <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.sign.up.tsv")
genes.down <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.sign.down.tsv")
allOE_genes <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.tsv")
#---------

# assign them new names
sign_genes_up <- genes.up$Symbol
sign_genes_down <- genes.down$Symbol
allgenes <- allOE_genes$Symbol

#perform enrichment analysis
ego.up <- enrichGO(gene = sign_genes_up, 
                   universe = allgenes,
                   keyType = "SYMBOL",  #ENSEMBL
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

ego.down <- enrichGO(gene = sign_genes_down, 
                     universe = allgenes,
                     keyType = "SYMBOL",  #ENSEMBL
                     OrgDb = org.Hs.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)


# Output results from GO analysis to a table
# take top 20 by gene count
top_go_up <- ego.up@result %>%
  arrange(desc(Count)) %>%
  head(20)

top_go_down <- ego.down@result %>%
  arrange(desc(Count)) %>%
  head(20)



# dotplot
p1 <- ggplot(top_go_up, aes(x = RichFactor, 
                         y = reorder(Description, RichFactor), 
                         size = Count, 
                         color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +#alpha = 0.8
  scale_color_gradientn(
    name = "log(BH adj.p)",
    colors = c("#fee8c8",  "#e34a33","#99000D"),# light→medium→dark
    limits = c(0.5, 5),                             # <-- define scale range
    breaks = c(0.5,2.5, 5),
    oob = scales::squish                           # clamp values outside limits
  )+
  theme_bw(base_size = 14) +
  labs(
    x = "Rich Factor",
    y = "GO term",
    color = expression(-log[10]("BH adj.p")),
    size = "Gene Count"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "italic")
  )
p1

p2 <- ggplot(top_go_down, aes(x = RichFactor, 
                            y = reorder(Description, RichFactor), 
                            size = Count, 
                            color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +#alpha = 0.8
  scale_color_gradientn(
    name = "log(BH adj.p)",
    colors = c("#a6cee3","#1f78b4", "darkblue"),
    limits = c(0.5, 5),                             # <-- define scale range
    breaks = c(0.5,2.5, 5),
    oob = scales::squish                           # clamp values outside limits
  )+
  theme_bw(base_size = 14) +
  labs(
    x = "Rich Factor",
    y = "GO term",
    color = expression(-log[10]("BH adj.p")),
    size = "Gene Count"
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, face = "italic")
  )
p2


## CHORD PLOTS -----------------------------------------------------------------
## Load datasets and perform the analysis independently for both timepoints
#---------
sigOE <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_6h_annotated.sign.tsv")
sigOE <- read.delim("~/Desktop/Analysis H.pylori/DESeq2_WT_vs_Mock_24h_annotated.sign.tsv")
#---------

allgenes <- allOE_genes$Symbol 
sigOE_genes <- as.character(sigOE$Symbol)

#perform enrichment analysis
ego <- enrichGO(gene = sigOE_genes, 
                universe = allgenes,
                keyType = "SYMBOL",  #ENSEMBL
                OrgDb = org.Hs.eg.db,  
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)



# Extract the enriched GO terms
go_gene_links <- as.data.frame(ego)
go_gene_links <- go_gene_links[c(1,2,9,11,12)] 

# Filter the top 5 GO terms by the smallest adjusted p-value
top_GO_terms <- go_gene_links %>%
  group_by(Description) %>%
  summarise(min_p = min(p.adjust)) %>%  # Find the minimum p-value for each term
  arrange(min_p) %>%  # Sort GO terms by p-value
  slice_head(n = 5) %>%  # Select the top 5
  pull(Description)  # Extract the GO term names

#--->>>>
#filter by counts
top_GO_terms <- go_gene_links %>%
  group_by(Description) %>%
  summarise(total_count = sum(Count)) %>%  # or max(Count), depending on how your df is structured
  arrange(desc(total_count)) %>%           # sort by highest Count
  slice_head(n = 5) %>%                    # take top 5
  pull(Description)
#<----

# Filter the original GO results for the top 5 GO terms
go_gene_links <- go_gene_links %>%
  filter(Description %in% top_GO_terms) %>%
  group_by(Description) %>%
  ungroup()

library(tidyr)
# Step 1: Split the concatenated gene IDs into individual rows
go_gene_links <- go_gene_links %>%
  separate_rows(geneID, sep = "/")  # Split by '/' to separate individual gene IDs

# To reduce the number of genes per term (example: take top 10 for each term)
go_gene_links <- go_gene_links %>%
  group_by(Description) %>%
  slice_head(n = 15) %>%  # Choose the top 15 genes for each term
  ungroup()

# Check the filtered data
head(go_gene_links)

#go_gene_links <- go_gene_links[-c(1)]#use pval
#go_gene_links <- go_gene_links[c(1,3,2)]
description <- go_gene_links # keep the table
go_gene_links <- go_gene_links[-c(2,3)] # use counts and GO terms
#go_gene_links <- go_gene_links[-c(1,3)] # use counts and description

# Step 6: Create the chord diagram with the filtered data
chordDiagram(go_gene_links, 
             transparency = 0.5,
             annotationTrack = c("name", "grid"))  # Pre-allocate space to avoid overlapping issues


#use tilted labels!
data <- go_gene_links

#create a chord diagram but without labeling 
chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 1)

#add the labels and axis
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 2.5, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.9)
  
  #print axis 
  #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, 
  #            sector.index = sector.name, track.index = 2)
  
}, bg.border = NA)



##### WITH GENE colors!!! ######
# join log2.FC into the GO-gene links
genes <- sigOE[c(2,5)]
colnames(genes) <- c("geneID","logFC")

go_gene_links <- go_gene_links %>%
  left_join(genes, by = "geneID") %>%   # change "geneID" if column is called "gene.name"
  mutate(direction = ifelse(logFC > 0, "up", "down"))


#define colors
link_colors <- ifelse(go_gene_links$direction == "up", "red", "steelblue")

# Create a named vector for sector colors
# GO terms = grey, genes = red/blue
sector_colors <- c(
  setNames(rep("grey70", length(unique(go_gene_links$Description))),
           unique(go_gene_links$Description)),
  setNames(ifelse(go_gene_links$direction == "up", "red", "blue"),
           go_gene_links$geneID)
)

# Step 1: basic chord diagram without labels
chordDiagram(go_gene_links,
             transparency = 0.5,
             annotationTrack = "grid",
             preAllocateTracks = 1,
             grid.col = sector_colors)  # add colors

# Step 2: add tilted labels like before
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  
  circos.text(mean(xlim), ylim[1] + 2.5, sector.name,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5),
              cex = 1)  # <-- increase this value to make font bigger
}, bg.border = NA)



##########################
# Prepare the description column
description <- description[c(1,2)]
description <-  description %>% distinct()

go_gene_links <- go_gene_links %>%
  left_join(description, by = "ID")

library(stringr)
go_gene_links <- go_gene_links %>%
  mutate(
    GO_label = str_wrap(Description, width = 15)
  )

# Create sector colors
# Draw chord diagram without sector.width
par(mar = c(1, 1, 1, 1))  # bottom, left, top, right margins
chordDiagram(go_gene_links[, c("GO_label","geneID")],
             transparency = 0.5,
             annotationTrack = "grid",
             preAllocateTracks = 1,
             grid.col = sector_colors)  # <-- remove sector.width

# Add tilted labels
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), 
              ylim[1] + 2.5, 
              sector.name,
              facing = "clockwise", 
              niceFacing = TRUE,
              adj = c(0, 0.5),
              cex = 0.8)
}, bg.border = NA)


