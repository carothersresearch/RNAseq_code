
##DESEQ

#Set working directory and load DESeq libraries, after installing them the first time. The path to the libraries is LOCAL, so you'll need to change it to point to your local installation. You can just use the checkboxes on the right to load the packages instead.
setwd("Desktop/")
library("DESeq", lib.loc="C:/Program Files/R/R-3.0.2/library")
library("DESeq2", lib.loc="C:/Program Files/R/R-3.0.2/library")

#Input tab-delimited text file, with first column gene names, all other columns HTSeq counts. One row of headers needs to be included.
countsTable <- read.delim("deseq_input.txt",header=TRUE, row.names=1)

#Define columns into groups
ODColData <- data.frame(MAIN=factor(c("EEU_1", "EEU_2", "EEU_3", "EEU_5", "EEU_7", "H2_1", "H2_2", "H2_4", "H2_6", "SWT_1", "SWT_2", "SWT_3", "SWT_4", "SWA_1", "SWA_2", "SWA_3", "SWA_4", "MB_1", "MB_3", "MB_4")),
                        row.names=colnames(countsTable))
#Make "other" groups
ODColData <- data.frame(MAIN=factor(c("EEU", "EEU", "EEU", "EEU", "EEU", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "MB", "MB", "MB")),
                        row.names=colnames(countsTable))
#Just run this
OD_dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=ODColData, design = ~ MAIN )
OD_dds_results <- DESeq(OD_dds, betaPrior=FALSE)

#List results; If your comparison of interest doesn't show up, then re-order the columns in your input file and run again. For some reason it doesn't give you the option to compare every pairwise sample group.
resultsNames(OD_dds_results)

#sub in the "name" value with a value in the results list
result <- results(OD_dds_results, name="MAIN_Other_vs_EEU")

#Enter filename to write to
write.table(result, file="MAIN_Other_vs_EEU", sep="\t")

#use Contrast to make specific comparisons, MAIN defines columns and EEU and Other are the comparisons to make
result_with_contrast <- results(OD_dds_results, contrast=c("MAIN","EEU","Other"))
#Write table
write.table(result_with_contrast, file="output.txt", sep="\t")


