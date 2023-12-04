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
alphadiversity <- physeqtax %>%                                                              
  plot_richness(
    x = "Treatment",                                                
    measures = c("Shannon", "Simpson")) +                           
  geom_violin(aes(fill = Treatment), show.legend = FALSE) +          
  geom_boxplot(width=0.1) +                                          
  theme_classic()+                                                   
  xlab(NULL)+                                                        
  theme(axis.text.y.left = element_text(size = 15),                  
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           
        axis.title.y = element_text(size = 15))+                     
  theme(strip.text = element_text(face = "bold", size = 20))+        
  scale_fill_manual(values=c("lightgoldenrod","darkorange", "darkred"))+   
  ggtitle("Alpha Diversity (Bacteria)") +                                       
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
  ordination = all_pcoa)+                                               
  geom_point(aes(fill = Treatment, shape = SampleType), size = 3) +     
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values=c("lightgoldenrod","darkorange", "darkred")) +
  theme_classic() +                                                     
  theme(                             
    legend.text = element_text(size = 20),                             
    legend.title = element_blank())+                                  
  theme(axis.text.y.left = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))        

#NMDS
nmds<- ordinate(physeqtax, "NMDS", "bray")
nmds_plot = plot_ordination(physeqtax, nmds, type="taxa", color="Phylum", title="taxa")+
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values=c("lightgoldenrod","darkorange", "darkred")) + 
  theme(                             
    legend.text = element_text(size = 15),                       
    legend.title = element_blank())+                   
  theme(axis.text.y.left = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))     
print(nmds_plot)
all_nmds <- nmds_plot + facet_wrap(~Phylum, 4)

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
  select(Phylum, Treatment, Abundance, SampleType) %>%
  group_by(Phylum, Treatment, SampleType) %>%
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
# Create a factor for SampleType to control the display of x-axis labels
avg.abundance$SampleType <- factor(avg.abundance$SampleType, levels = unique(avg.abundance$SampleType))

phylum_avgabundance_sampletype <- ggplot(avg.abundance, aes(x = Treatment, y = avg_abundance, fill = Phylum)) +
  geom_col(position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL) +
  facet_wrap(~SampleType, ncol = 1) + 
  theme(strip.text.x = element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10), 
    strip.placement = "outside",  
    title = element_text(size = 10)
  )

#Calcultae proportion of each phyla
phylumabundance <- as.data.frame(phylumabundance)
genus_colors <- c(
  "#5F7FC7", "orange","#DA5724", "#508578",
  "#AD6F3B", "#673770", "#652926", "#C84248",  "#D1A33D")

#proportion of the community "Acidobacteriota"
Acidobacteriota <- subset(phylumabundance, Phylum =="Acidobacteriota")
Acidobacteriota_avgabundance <- ggplot(Acidobacteriota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Acidobacteriota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  
  facet_wrap(~ SampleType, ncol = 1)

print(Acidobacteriota_avgabundance)

#proportion of the community "Actinobacteriota"
Actinobacteriota <- subset(phylumabundance, Phylum =="Actinobacteriota")
Actinobacteriota_avgabundance <- ggplot(Actinobacteriota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Actinobacteriota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  
  facet_wrap(~ SampleType, ncol = 1)

print(Actinobacteriota_avgabundance)

#proportion of the community "Bacteroidota"
bacteroidota <- subset(phylumabundance, Phylum =="Bacteroidota")
bacteroidota_avgabundance <- ggplot(bacteroidota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Bacteroidota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  
  facet_wrap(~ SampleType, ncol = 1)

print(bacteroidota_avgabundance)
#proportion of the community "Chloroflexi" 
Chloroflexi <- subset(phylumabundance, Phylum =="Chloroflexi")
Chloroflexi_avgabundance <- ggplot(Chloroflexi) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Chloroflexi") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  
  facet_wrap(~ SampleType, ncol = 1)

print(Chloroflexi_avgabundance)

#proportion of the community "Gemmatimonadota"
Gemmatimonadota <- subset(phylumabundance, Phylum =="Gemmatimonadota")
Gemmatimonadota_avgabundance <- ggplot(Gemmatimonadota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Gemmatimonadota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) + 
  facet_wrap(~ SampleType, ncol = 1)

print(Gemmatimonadota_avgabundance)

#proportion of the community "Planctomycetota"
Planctomycetota <- subset(phylumabundance, Phylum =="Planctomycetota")
Planctomycetota_avgabundance <- ggplot(Planctomycetota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Planctomycetota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  
  facet_wrap(~ SampleType, ncol = 1)
print(Planctomycetota_avgabundance)

#proportion of the community "Verrucomicrobiota"
Verrucomicrobiota <- subset(phylumabundance, Phylum =="Verrucomicrobiota")
Verrucomicrobiota_avgabundance <- ggplot(Verrucomicrobiota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Verrucomicrobiota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) + 
  facet_wrap(~ SampleType, ncol = 1)

print(Verrucomicrobiota_avgabundance)

#heatmap
#Phyla
s <- plot_heatmap(physeqtax,  method = "NMDS", distance = "bray", 
                  taxa.label = "Phylum", taxa.order = "Phylum", 
                  low="darkblue", high="red", na.value="white")+ 
  theme(axis.text.y.left = element_text(size = 5),
        axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 5),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20)))

s 
#Heatmap for a subset of samples
new_palette <- colorRampPalette(c("darkblue", "white", "red"))

plot_taxa_heatmap(physeqtax, subset.top = 50,
                                   VariableA = "Treatment",
                                   heatcolors = new_palette(100),
                                   transformation = "log10")




