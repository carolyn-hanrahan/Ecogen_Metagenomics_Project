## Ecological Genomics Group Project Notebook

### November - December 2023

### Alex, Noah, Lucy, Maryam and Carolyn

### Group Work Day #1: November 27th, 2023

-   Remove columns with count numbers of 1 or 2 in the `its_sequence_data.xlsx` file and `16S_sequence_data.xlsx` file
-   List of graphics and analyses we are interested in implementing:
    -   Pie charts
        -   comparing microbes across root locations
    -   Baysien clustering
        -   comparing fertilizer treatment
    -   Heatmap
        -   correlating taxa --> phenotypic data
        -   indicator species analysis
    -   WGCNA
    -   Functional gene analysis of specific/dominant taxa
        -   Piecrust package in R
- Imported code into R and began working on manipulating the data and following a tutorial using the `phyloseq` library. Link to tutorial here: 

https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#step-by-step 

### Group Work Day #2: November 29th, 2023

- Maryam successfully implemented phyloseq program with our data!
- Examined mctoolsr scripts on Github: https://github.com/leffj/mctoolsr 
- Copied code into new script titled: `Script_Nov_29_2023.R`

- See `16S_bioconductor.R` file and `ITS_bioconductor.R` file. The corresponding excel sheets used for each respective analysis were: `16Sdata.xlsx` and `ITS_data.xlsx`. These files can be found in the data folders within the Ecogen_Metagenomics_Project. With this code, we successfully created some good bargraphs, heatmaps and plots of the various microbial taxa within our samples. 

### Group Work Day #3: December 4th, 2023 

- Added Maryam's R code `Shared_Meta16S.R` and `16Sdata.xlsx` file (which the R script is based off of)
- Examined rarefaction parameters + cutoffs. Created a plot to determine the cutoff sample size whereby the values even out. 
- Started with 90 samples and after rarefaction were left with 52 samples. 
- Created diversity graphs and barplots 
- Maintain statistical power and minimize bias 
- Altering rarefaction parameters. 
  - Cutoff of 180 to 10,000 
- Preliminary results: sampling location on the plant was more influential of a factor than treatment type (fertilizer)
- Next steps... 
  - Alter rarefaction metrics 
  - Pick the most abundant taxa and examine at the genus level... or identify phylums of bacteria that are most interesting
  - ANCOMBC package to look for statistical significance. This program will tell which taxa/ASVs are statistically different between the treatment groups [if time allows]
  - Use Shannon's diversity rather than other diversity metric 
  
Meet to work on writeup tomorrow: 1:30 PM! 
Wednesday: examine specific taxa 


