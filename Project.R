install.packages("htmltools")
library(htmltools)
library( DESeq2 )
library(ggplot2)

data <- read.table("D:/BioinfoProject/GSE184072_count_Data.txt", header = TRUE, row.names = 1)
head(data)
rownames(data)
colnames(data)
# Define the row names you want to delete
rows_to_delete <- c("1-Dec","1-Mar","1-Sep","10-Mar","10-Sep","11-Mar","11-Sep","12-Sep","14-Sep","2-Mar","2-Sep","3-Mar","3-Sep","4-Mar", "4-Sep","5-Mar","5-Sep", "6-Mar","6-Sep","7-Mar" , "7-Sep", "8-Mar", "8-Sep", "9-Mar","9-Sep" )

# Delete rows with the specified row names
data <- data[!(rownames(data) %in% rows_to_delete), ]

head(data)

colData<-read.csv("D:/BioinfoProject/sample_info.csv")
head (colData)
rownames(colData)
colnames(colData)



# making sure the row names in colData matches to column names in counts_data
all(colnames(data) %in% rownames(colData))

# are they in the same order?
all(colnames(data) == rownames(colData))

summary(data)
data <- round(data)
# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~ Tissue)

dds


# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$Tissue <- relevel(dds$Tissue, ref = "Normal")

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res

# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

png("D:/BioinfoProject/MAPlot.png", width = 800, height = 600)  
plotMA(res)
dev.off() 

png("D:/BioinfoProject/VolcanoPlot.png", width = 800, height = 600)  # Specify file name and dimensions
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()  # Save the plot
