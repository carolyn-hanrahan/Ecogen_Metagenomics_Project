setwd("/Users/maryamnouri-aiin/Desktop/CompostMetagenomics/16S/")

#Required packages
require(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(wesanderson)
library(RColorBrewer)
library("viridis")  
library(scales)
library(igraph)
library(ggridges)
library(agricolae)
library(patchwork)
library(qiime2R)
library(phyloseq)
library(BiocManager)
library(vegan)
library("ape")
library("DESeq2")
####################################Rarefaction and normalisation####################################
#Rawdata
data <- read.xlsx("16S_1.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE)
# Exclude the specified column "taxonomy
column_to_exclude <- c("taxonomy")
data <- data[, !colnames(data) %in% column_to_exclude]
data <- data[rowSums(data)>600,]
data <- data[,colSums(data)>19000]
OTU <- specnumber(t(data))
# Determine the rarefaction level
raremax <- min(colSums(data))
# Perform rarefaction
ACrare <- rarefy(t(data), raremax)
#Plot rarified data
plot(OTU, ACrare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# Customize the rarefaction curve
AC_Curve <- rarecurve(t(data), step = 20, sample = raremax, col = "blue", cex = 0.6)
# Perform relative rarefaction
AC_Data_rare <- rrarefy(t(data), raremax)
#.csv file to use for the rest of analysis
write.csv(t(AC_Data_rare), file = "AC_16S_Rfy.csv")

####################################https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#step-by-step####################################
####################################Phyloseq####################################
#Three tables are needed: OTU, Taxonomy, Samples
otu_mat<- as.data.frame(read_excel("16Sdata.xlsx", sheet = "OTU matrix"))
tax_mat<- read_excel("16Sdata.xlsx", sheet = "Taxonomy table")
samples_df <- as.data.frame(read_excel("16Sdata.xlsx", sheet = "Samples"))

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV #")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV #")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#Find the unique phylum and taxa
otumat <- data.frame(otu_mat)
taxmat <- data.frame(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq0 <- phyloseq(OTU, TAX)
physeq <- phyloseq(OTU, TAX, samples)
physeq

#Rarefy_even_depth: Resample an OTU table such that all samples have the same
ps.rarefied = rarefy_even_depth(physeq, rngseed=1, 
                                sample.size=0.5*min(sample_sums(physeq)), replace=F)

#Keep only specific phylum
physeqtax <- subset_taxa(ps.rarefied, Phylum %in% c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", 
                                                    "Chloroflexi", "Gemmatimonadota", "Planctomycetota", "Verrucomicrobiota"))
#Exclude classes
physeqtax <- subset_taxa(physeqtax, !(Class %in% c("BD7-11", "FFCH5909")))
# Exclude NA values for Genus and Species
physeqtax <- subset_taxa(physeqtax, !is.na(Genus) & !is.na(Species))
# Exclude Genus 
physeqtax <- subset_taxa(physeqtax, !(Genus %in% c("JG30-KF-CM45", "CPla-3 termite group", "67-14", "SRB2", "BIrii41", "env.OPS 17", 
                                                   "WWH38", "AKYG1722", "Amb-16S-1323", "A4b", "TRA3-20", "Unknown Family", "NA")))
# Exclude Species 
physeqtax <- subset_taxa(physeqtax, !(Species %in% c("Ellin516", "DEV008",  "SH-PL14", "JGI 0001001-H03", "BIyi10", "JCM 18997", 
                                                     "Ellin517", "NA")))
physeqtax


###Create violin plot of alpha diversity
alphadiversity_treatment_sampletype <- physeqtax %>%                                                            
  plot_richness(
    x = "Treatment",
    measures = c("Shannon")) +                      
  geom_violin(aes(fill = Treatment), show.legend = TRUE) +         
  geom_boxplot(width=0.1) +  
  facet_grid(PlantType~SampleType) +                               
  theme_classic()+                                                   
  xlab(NULL)+                                                       
  theme(axis.text.y.left = element_text(size = 15),                 
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           
        axis.title.y = element_text(size = 15))+                     
  theme(strip.text = element_text(face = "bold", size = 20))+        
  scale_fill_manual(values=c("lightgoldenrod","darkorange", "darkred"))+   
  ggtitle("") +                                        #add title
  theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))

#Making a PCoA plot

all_pcoa <- ordinate(
  physeq = physeqtax, 
  method = "PCoA", 
  distance = "jaccard"
)
#Plot
pca <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa) +
  geom_point(aes(fill = Treatment, shape = SampleType, color = PlantType), size = 25) +  
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("lightgoldenrod", "darkorange", "darkred")) +
  scale_color_manual(values = c("green", "red")) +  
  scale_size_manual(values = c(3, 4, 5)) +  
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.key.size = unit(2, "lines"))     


#Filter out eukaryotes and mitochondria
justbacteria <- physeqtax %>%
  subset_taxa(
    Kingdom == "Bacteria" 
  )
justbacteria

#Convert phyloseq object into a data frame with relative abundance counts
phylumabundance <- justbacteria %>%
  tax_glom(taxrank = "Genus") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  arrange(Phylum) 
head(phylumabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
avg.abundance <- phylumabundance %>%
  select(Phylum, Treatment, Abundance, PlantType, SampleType) %>%
  group_by(Phylum, Treatment, PlantType, SampleType) %>%
  summarize(
    avg_abundance = mean(Abundance), .groups = 'drop'
  ) %>%
  filter(avg_abundance > 0.0003)
head(avg.abundance)
unique(avg.abundance$Phylum)
dim(avg.abundance)

###Save phylum colors
phylum_colors <- c(
  "#5F7FC7", "orange","#DA5724", "#508578",
  "#AD6F3B", "#673770", "#652926", "#C84248",  "#D1A33D"
)
#heatmap
heatmap <- ggplot(avg.abundance, aes(x = Treatment, y = Phylum, fill = avg_abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgoldenrod", high = "darkred") +
  facet_grid(PlantType ~ SampleType) +  # Add facets for PlantType and SampleType
  theme_minimal() +
  labs(title = "Average Abundance of Phyla",
       x = "Treatment",
       y = "Phylum",
       fill = "Average Abundance")
# Create a factor for SampleType to control the display of x-axis labels
avg.abundance$SampleType <- factor(avg.abundance$SampleType, levels = unique(avg.abundance$SampleType))
avg.abundance$PlantType <- factor(avg.abundance$PlantType, levels = unique(avg.abundance$PlantType))

phylum_avgabundance_sampletype <- ggplot(avg.abundance, aes(x = Treatment, y = avg_abundance, fill = Phylum)) +
  geom_col(position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL) +
  facet_grid(PlantType ~ SampleType) +   
  theme(strip.text.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10),  # Adjust facet label size
    strip.placement = "outside",  # Place the facet labels outside the plot area
    title = element_text(size = 10)
  ) +
  labs(x = "Treatment", subtitle = "")  # Adding a subtitle for PlantType

phylum_avgabundance_sampletype

#Calcultae proportion of each phyla
phylumabundance <- as.data.frame(phylumabundance)
phylumabundance$SampleType <- factor(phylumabundance$SampleType, levels = unique(phylumabundance$SampleType))
phylumabundance$PlantType <- factor(phylumabundance$PlantType, levels = unique(phylumabundance$PlantType))
#e.g. proportion of the community "Acidobacteriota"
Acidobacteriota <- subset(phylumabundance, Phylum =="Acidobacteriota")
Acidobacteriota_avgabundance <- ggplot(Acidobacteriota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL) +
  facet_grid(PlantType ~ SampleType) +   
  theme(strip.text.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10),  # Adjust facet label size
    strip.placement = "outside",  # Place the facet labels outside the plot area
    title = element_text(size = 10)
  ) +
  labs(x = "Treatment", subtitle = "")  # Adding a subtitle for PlantType


print(Acidobacteriota_avgabundance)

#heatmap
#Phyla
s <- plot_heatmap(physeqtax,  method = "NMDS", distance = "bray", 
             taxa.label = "Phylum", taxa.order = "Phylum", 
             low="lightgoldenrod", high="darkred", na.value="white")
s 
#Heatmap for a subset of samples
new_palette <- colorRampPalette(c("lightgoldenrod", "white", "darkred"))


heatmap_top50 <- plot_taxa_heatmap(physeqtax,
                                 subset.top = 50,
                                 VariableA = "Treatment",
                                 heatcolors = new_palette(100),
                                 transformation = "log10")
heatmap_top50




