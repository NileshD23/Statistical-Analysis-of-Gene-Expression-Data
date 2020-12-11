# A. install DEseq 
# -- Warning, this will take ~5 mins to install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
# B. Add Code to read in RNA-Seq data

###################
#these lines need to be updated with the appropriate filenames!!
###################
BCRS <- read.table("BC_RNAseq.txt",header=T,row.names=1,sep="\t",quote="")
BCRS_stat <- read.table("BC_RNAseq_Status.txt",header=T,sep="\t",quote="")
###################


rownames(BCRS_stat) <- BCRS_stat$BCID
# C. Convert Status table into a factor
BCRS_stat$ERStat <- factor(BCRS_stat$ERStat)
cds <- DESeqDataSetFromMatrix(countData = BCRS,
                              colData = BCRS_stat,
                              design = ~ ERStat)
# D. Normalize for Size Effects:
cds <- DESeq(cds)

# E. Calculate Differential Expression
res <- results(cds,contrast = c('ERStat','0','1'))

# F. Extract top 100 genes by their adusted p-value using FDR
res <- res[order(res$padj), 1:100]
print(res)

print(FDR(data = res))