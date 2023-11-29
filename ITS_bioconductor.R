setwd("C:/Users/Carolyn/Desktop/Ecogen Group Project/Ecogen_Metagenomics_Project/ITS_data") 

library(tidyverse)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(RColorBrewer)
library("viridis")  
library(scales)
library(agricolae)
require(readxl)
library(phyloseq)
library(patchwork)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(MicrobiotaProcess)
library(BiocManager)
BiocManager::version()
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install(version = "devel")
BiocManager::install("vcfR", force = TRUE)
BiocManager::valid()
source("https://bioconductor.org/biocLite.R")
BiocManager::install("phyloseq", force = TRUE)
library("phyloseq")
source("http://bioconductor.org/biocLite.R")
biocLite("multtest", type="source")
##################
################
devtools::install_github("joey711/phyloseq")
#################
###############

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

install(
  pkgs = character(),
  site_repository = character(),
  update = TRUE,
  ask = TRUE,
  checkBuilt = FALSE,
  force = FALSE,
  version = BiocManager::version()
)
BiocManager::install("phyloseq", type="source")

bacteria <- read_excel("ITS_data.xlsx")
otu_mat<- read_excel("ITS_data.xlsx", sheet = "OTU Matrix")
tax_mat<- read_excel("ITS_data.xlsx", sheet = "Taxonomy Table")
samples_df <- read_excel("ITS_data.xlsx", sheet = "Samples")

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV #") 

tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV #")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample #") 

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
taxmat <- data.frame(tax_mat)
is.atomic(taxmat$Class)
unique(taxmat$Class)
unique(taxmat$Phylum)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)


carbom <- phyloseq(OTU, TAX, samples)
carbom
sample_names(carbom)

rank_names(carbom)
sample_variables(carbom)


carbom <- subset_taxa(carbom, Phylum %in% c("Bacteroidota", "Acetothermia", "Acidobacteriota", 
                                            "Aquificota", "Armatimonadota", "Chloroflexi"))
carbom <- subset_taxa(carbom, !(Class %in% c("Abditibacteria", "Vicinamibacteria")))
carbom


total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)


plot_bar(carbom, fill = "Phylum")

plot_bar(carbom, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

carbom_fraction <- merge_samples(carbom, "Treatment")
plot_bar(carbom_fraction, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


carbom_bact <- subset_taxa(carbom, Phylum %in% c("Bacteroidota"))
plot_bar(carbom_bact, x="Genus", fill = "Genus", facet_grid = Treatment~SampleType) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

plot_heatmap(carbom, method = "NMDS", distance = "bray")

carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.20) > 0, TRUE)
carbom_abund
otu_table(carbom_abund)[1:8, 1:5]
plot_heatmap(carbom_abund, method = "NMDS", distance = "bray")

plot_heatmap(carbom_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class", 
             trans=NULL, low="beige", high="red", na.value="beige")

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

plot_heatmap(carbom_bact, method = "NMDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             low="beige", high="red", na.value="beige")

plot_richness(carbom, measures=c("Chao1", "Shannon"))

plot_richness(carbom, measures=c("Chao1", "Shannon"), x="Treatment", color="SampleType")

carbom.ord <- ordinate(carbom, "NMDS", "bray")

plot_ordination(carbom, carbom.ord, type="taxa", color="Class", shape= "Phylum", 
                title="ASV")


plot_ordination(carbom, carbom.ord, type="split", color="Class", 
                shape="Treatment", title="biplot", label = "SampleType") +  
  geom_point(size=3)

plot_net(carbom, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="Class", point_label="Genus")

plot_net(carbom_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Class", point_label="Genus") 
