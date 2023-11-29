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
