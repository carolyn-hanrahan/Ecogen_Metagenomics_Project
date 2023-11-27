# Notes from November 27, 2023
# Carolyn Hanrahan 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Set your working directory
setwd("C:/Users/Carolyn/Desktop/Ecogen Group Project/Ecogen_Metagenomics_Project/ITS_data") 


#install.packages("readxl")
library(readxl)
library(dplyr)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq") # Install phyloseq


install.packages(c("RColorBrewer", "patchwork")) #install patchwork to chart publication-quality plots and readr to read rectangular datasets.

library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Specify the file path for ITS data 
ITS_file_path <- "C:/Users/Carolyn/Desktop/Ecogen Group Project/Ecogen_Metagenomics_Project/ITS_data/ITS_Sequence_Data_clean.xlsx"

# Read the Excel file & Import ITS data
its_data <- read_excel(ITS_file_path)

head(its_data)


#Specify the file path for 16S data
sixteen_s_file_path <- "C:/Users/Carolyn/Desktop/Ecogen Group Project/Ecogen_Metagenomics_Project/16S_data/16S_Sequence_Data_clean.xlsx"

# Read the Excel file and import 16S data: 
sixteen_s_data <- read_excel(sixteen_s_file_path)

head(sixteen_s_data)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Practice Code

# converting ITS excel file into a dataframe: 

its_df <- as.data.frame(its_data)
head(its_df)

# subsetting dataframe to examine sample numbers and corresponding treatments:
samples_and_treatment <- its_df[ ,1:2]

# converting 16S excel file into a dataframe: 

sixteen_s_df <- as.data.frame(sixteen_s_data)

head(sixteen_s_df)





