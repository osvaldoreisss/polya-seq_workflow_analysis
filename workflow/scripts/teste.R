library("DESeq2")                                                            
library("ggplot2")
library("RColorBrewer")
library("pheatmap")

samples_table <- read.csv("../../config/samples.tsv", header = TRUE, sep=",")


files_three <- paste(samples_table$sample, "_three_out.txt", sep="")
files_five <- paste(samples_table$sample, "_five_out.txt", sep="")
files_cds <- paste(samples_table$sample, "_cds_out.txt", sep="")

counts_table_three = matrix(, nrow = 7072, ncol = 0)
for(f in files_three){
    tmp = read.table(f, header=F)
    colnames(tmp) = c("gene", f)
    counts_table_three = cbind(counts_table_three, tmp[order(tmp[,1]), , drop=F])
}
rownames(counts_table_three) <- counts_table_three$gene
counts_table_three = counts_table_three[, colnames(counts_table_three) != "gene"]

colData <- data.frame(row.names=colnames(counts_table_three), condition=as.factor(samples_table$condition))
dds_three <- DESeqDataSetFromMatrix(countData=counts_table_three, colData=colData, design = ~ condition )

dds_three <- DESeq(dds_three)
res_three <- results(dds_three, contrast=c("condition","YPD_diauxic_yeast_cells","YPD_log_phase_yeast_cells"))
res_threeOrdered <- res_three[order(res_three$pvalue),]
head(res_threeOrdered)


counts_table_five = matrix(, nrow = 7072, ncol = 0)
for(f in files_five){
    tmp = read.table(f, header=F)
    colnames(tmp) = c("gene", f)
    counts_table_five = cbind(counts_table_five, tmp[order(tmp[,1]), , drop=F])
}
rownames(counts_table_five) <- counts_table_five$gene
counts_table_five = counts_table_five[, colnames(counts_table_five) != "gene"]

colData <- data.frame(row.names=colnames(counts_table_five), condition=as.factor(samples_table$condition))
dds_five <- DESeqDataSetFromMatrix(countData=counts_table_five, colData=colData, design = ~ condition )

dds_five <- DESeq(dds_five)
res_five <- results(dds_five, contrast=c("condition","YPD_diauxic_yeast_cells","YPD_log_phase_yeast_cells"))
res_fiveOrdered <- res_five[order(res_five$pvalue),]
head(res_fiveOrdered)

counts_table_cds = matrix(, nrow = 7072, ncol = 0)
for(f in files_cds){
    tmp = read.table(f, header=F)
    colnames(tmp) = c("gene", f)
    counts_table_cds = cbind(counts_table_cds, tmp[order(tmp[,1]), , drop=F])
}
rownames(counts_table_cds) <- counts_table_cds$gene
counts_table_cds = counts_table_cds[, colnames(counts_table_cds) != "gene"]


colData <- data.frame(row.names=colnames(counts_table_cds), condition=as.factor(samples_table$condition))
dds_cds <- DESeqDataSetFromMatrix(countData=counts_table_cds, colData=colData, design = ~ condition )

dds_cds <- DESeq(dds_cds)
res_cds <- results(dds_cds, contrast=c("condition","YPD_diauxic_yeast_cells","YPD_log_phase_yeast_cells"))
res_cdsOrdered <- res_cds[order(res_cds$pvalue),]
head(res_cdsOrdered)
