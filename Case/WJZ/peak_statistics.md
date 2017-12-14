m6A peak statistics
=======
### including:
1. venn diagram to plot peak overlapping across samples
2. calculating number of peaks per transcript
3. calculating how many m6A sites were rescued
4. compare m6A sites between samples
5. peak scatter plot (mimicing **_Nature Letters_** m6A modulates haematopoietic stem and progenitor cell specification)

1 peak venn diagram
=====
```R
ctrl_1 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_ctrl_1','con_peak.bed')
ctrl_2 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_ctrl_2','con_peak.bed')
ctrl_3 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_ctrl_3','con_peak.bed')

kd_1 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_KD_1','con_peak.bed')
kd_2 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_KD_2','con_peak.bed')
kd_3 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','IP_KD_3','con_peak.bed')

res_1 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','res_rep_1','con_peak.bed')
res_2 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','res_rep_2','con_peak.bed')
res_3 <- file.path('D:','RWD','WJZ','m6a-seq','RMT_peaks','res_rep_3','con_peak.bed')




library(ChIPseeker)
library(grid)
library(VennDiagram)

a <- read.table(ctrl_1, sep = '\t')
a <- merge(a, read.table(ctrl_2, sep = '\t'), by = 4)
a <- merge(a, read.table(ctrl_3, sep = '\t'), by.x = 1, by.y = 4)
a <- a[!duplicated(a$V4),]
a <- a[,c(1,2)]

b <- read.table(kd_2, sep = '\t')
b <- merge(b, read.table(kd_3, sep = '\t'), by = 4)
b <- b[!duplicated(b$V4),]
b <- b[,c(1,2)]

c <- read.table(res_2, sep = '\t')
c <- merge(c, read.table(res_3, sep = '\t'), by = 4)
c <- c[!duplicated(c$V4),]
c <- c[,c(1,2)]

abc <- merge(a,b, by = 1)
abc <- merge(abc,c, by = 1)

ab <- merge(a, b, by = 1)
ac <- merge(a, c, by = 1)
bc <- merge(b, c, by = 1)

grid.newpage()
draw.triple.venn(area1 = 3289, area2 = 8822, area3 = 3527, n12 = 3240, n23 = 3446, n13 = 2597, 
                 n123 = 2581, category = c("replicate 1", "replicate 2", "replicate 3"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"), scaled = T, print.mode = c('raw','percent'), 
                 #cex = 0, cat.cex = 0
                 )


# re-enable the scaling funciton
# overrideTriple <- 1
grid.newpage()
draw.triple.venn(area1 = 6520, area2 = 6864, area3 = 6511, n12 = 4297, n23 = 3996, n13 = 4299, 
                 n123 = 3544, category = c("replicate 1", "replicate 2", "replicate 3"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"), scaled = T)


vennplot.peakfile(c(kd_2,kd_3))
grid.newpage()
# draw.triple.venn(area1 = 5871, area2 = 14161, area3 = 19591, n12 = 4700, n23 = 11901, n13 = 5026, 
#                  n123 = 4615, category = c("replicate 1", "replicate 2", "replicate 3"), lty = "blank", 
#                  fill = c("skyblue", "pink1", "mediumorchid"))
draw.pairwise.venn(area1 = 14161, area2 = 19591, cross.area = 11901,  
                 category = c("replicate 2", "replicate 3"), lty = "blank", 
                 fill = c("skyblue", "pink1"))


vennplot.peakfile(c(res_2,res_3))
grid.newpage()
# draw.triple.venn(area1 = 13724, area2 = 11265, area3 = 5230, n12 = 8359, n23 = 4280, n13 = 4048, 
#                  n123 = 3916, category = c("replicate 1", "replicate 2", "replicate 3"), lty = "blank", 
#                  fill = c("skyblue", "pink1", "mediumorchid"))
draw.pairwise.venn(area1 = 5230, area2 = 11265, cross.area = 4280,  
                   category = c("replicate 2", "replicate 3"), lty = "blank", 
                   fill = c("skyblue", "pink1"))


vennplot.peakfile(c(ctrl_1,ctrl_2, kd_2, kd_3, res_2, res_3))
draw.triple.venn(area1 = 6520, area2 = 6864, area3 = 6511, n12 = 4297, n23 = 3996, n13 = 4299, 
                 n123 = 3544, category = c("replicate 1", "replicate 2", "replicate 3"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"), scaled = T)
```

2 peak numbers per transcript
======
```R
peak1 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_ctrl_1\\con_peak.xls', header = T,
                   sep = '\t')
table1 <- as.data.frame(table(peak1$name))
plot(table(table1$Freq))

count1 <- sum(table1$Freq)/nrow(table1)

peak2 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_ctrl_2\\con_peak.xls', header = T,
                    sep = '\t')
table2 <- as.data.frame(table(peak2$name))
plot(table(table2$Freq))

count2 <- sum(table2$Freq)/nrow(table2)


peak3 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_ctrl_3\\con_peak.xls', header = T,
                    sep = '\t')
table3 <- as.data.frame(table(peak3$name))
plot(table(table3$Freq))

count3 <- sum(table3$Freq)/nrow(table3)


peak4 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_KD_1\\con_peak.xls', header = T,
                    sep = '\t')
table4 <- as.data.frame(table(peak4$name))
plot(table(table4$Freq))

count4 <- sum(table4$Freq)/nrow(table4)


peak5 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_KD_2\\con_peak.xls', header = T,
                    sep = '\t')
table5 <- as.data.frame(table(peak5$name))
plot(table(table5$Freq))

count5 <- sum(table5$Freq)/nrow(table5)


peak6 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_ctrl_kd_Result\\IP_KD_3\\con_peak.xls', header = T,
                    sep = '\t')
table6 <- as.data.frame(table(peak6$name))
plot(table(table6$Freq))

count6 <- sum(table6$Freq)/nrow(table6)


peak7 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_rescue_Result\\res_rep_1\\con_peak.xls', header = T,
                    sep = '\t')
table7 <- as.data.frame(table(peak7$name))
plot(table(table7$Freq))

count7 <- sum(table7$Freq)/nrow(table7)


peak8 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_rescue_Result\\res_rep_2\\con_peak.xls', header = T,
                    sep = '\t')
table8 <- as.data.frame(table(peak8$name))
plot(table(table8$Freq))

count8 <- sum(table8$Freq)/nrow(table8)


peak9 <- read.table('D:\\RWD\\WJZ\\exomePeak\\RMT_rescue_Result\\res_rep_3\\con_peak.xls', header = T,
                    sep = '\t')
table9 <- as.data.frame(table(peak9$name))
plot(table(table9$Freq))

count9 <- sum(table9$Freq)/nrow(table9)


count <- data.frame(count = c(count1, count2, count3, count4, count5, count6, count7, count8, count9))

plot(count$count)
mean(count$count)
```

3 peak rescue percentage
===
```R
load('WJZ/starburst_volcano_kd_ctrl.Rdata')
volcano_kd <- volcano_normal
rm(volcano_normal)
load('WJZ/starburst_volcano_res_kd.Rdata')
volcano_res <- volcano_normal
rm(volcano_normal)

volcano_all <- merge(volcano_kd, volcano_res, by = 'gene_id')
rm(volcano_kd, volcano_res)

write.csv(volcano_all, 'volcano_all.csv')

volcano_all <- volcano_all[,c(5,1,6,7,2,4,8,9,12,17,21,23,27,28,31,36)]
colnames(volcano_all) <- c('gene_symbol', 'gene_id', 'entrez_id', 'gene_biotype', 
                           'logFC.kd.m6a', 'FDR.kd.m6a', 'logFC.kd.mrna', 'logCPM.kd.mrna', 'FDR.kd.mrna', 'starburst.status.kd', 
                           'logFC.res.m6a', 'FDR.res.m6a', 'logFC.res.mrna', 'logCPM.res.mrna', 'FDR.res.mrna', 'starburst.status.res')


keep <- volcano_all$FDR.kd.m6a<0.05&volcano_all$FDR.kd.mrna<0.05&volcano_all$FDR.res.m6a<0.05&volcano_all$FDR.res.mrna<0.05
table(keep)
volcano_sig <- volcano_all[keep,]

FC_cut <- 0.5
keep2 <- abs(volcano_sig$logFC.kd.m6a)>FC_cut&abs(volcano_sig$logFC.kd.mrna)>FC_cut&abs(volcano_sig$logFC.res.m6a)>FC_cut&abs(volcano_sig$logFC.res.mrna)>FC_cut
table(keep2)
volcano_sig2 <- volcano_sig[keep2,]

keep_tar <- volcano_sig$logFC.kd.m6a > 0 
volcano_sig_tar <- volcano_sig[keep_tar,]

keep2_tar <- volcano_sig2$logFC.kd.m6a > 0 
volcano_sig2_tar <- volcano_sig2[keep2_tar,]

write.csv(volcano_sig2_tar, 'rescued gene.csv')

# num expected
expected <- (volcano_sig$logFC.kd.m6a*volcano_sig$logFC.res.m6a<0)&(volcano_sig$logFC.kd.mrna*volcano_sig$logFC.res.mrna<0)
table(expected)

expected2 <- (volcano_sig2$logFC.kd.m6a*volcano_sig2$logFC.res.m6a<0)&(volcano_sig2$logFC.kd.mrna*volcano_sig2$logFC.res.mrna<0)
table(expected2)

expected_tar <- (volcano_sig_tar$logFC.kd.m6a*volcano_sig_tar$logFC.res.m6a<0)&(volcano_sig_tar$logFC.kd.mrna*volcano_sig_tar$logFC.res.mrna<0)
table(expected_tar)

expected2_tar <- (volcano_sig2_tar$logFC.kd.m6a*volcano_sig2_tar$logFC.res.m6a<0)&(volcano_sig2_tar$logFC.kd.mrna*volcano_sig2_tar$logFC.res.mrna<0)
table(expected2_tar)

write.csv(volcano_sig2_tar[expected2_tar,],'rescued_gene.csv')

library(ggplot2)


library(clusterProfiler)

# data(geneList, package="DOSE")
# gene <- names(geneList)[abs(geneList) > 2]

# gene.df <- rescued[,c(6,1,5)]
# colnames(gene.df) <- c('ENTREZID', 'ENSEMBL', 'SYMBOL')

gene <- volcano_sig[expected,]$entrez_id

library(org.Hs.eg.db)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

barplot(ggo, drop=TRUE, showCategory=12)

ego <- enrichGO(gene          = gene,
                #universe      = volcano_all$entrez_id,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

barplot(ego, showCategory=15)
dotplot(ego,  showCategory = 15)
enrichMap(ego)


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
browseKEGG(kk, 'hsa04621')
```

4 compare m6A peaks
===
```R
kd_vs_ctrl <- read.table('D:\\RWD\\WJZ\\kd_vs_ctrl_exomePeak_out\\exomePeak_output\\con_sig_diff_peak.xls', 
                         header = T, stringsAsFactors = F, sep = '\t')
kd_vs_ctrl <- kd_vs_ctrl[,c('name', 'diff.log2.fc', 'diff.lg.p', 'diff.lg.fdr')]

oe_vs_ctrl <- read.table('D:\\RWD\\WJZ\\oe_vs_ctrl_exomePeak_out\\exomePeak_output\\con_sig_diff_peak.xls', 
                         header = T, stringsAsFactors = F, sep = '\t')
oe_vs_ctrl <- oe_vs_ctrl[,c('name', 'diff.log2.fc', 'diff.lg.p', 'diff.lg.fdr')]

oe_vs_kd <- read.table('D:\\RWD\\WJZ\\oe_vs_kd_exomePeak_out\\exomePeak_output\\con_sig_diff_peak.xls', 
                         header = T, stringsAsFactors = F, sep = '\t')
oe_vs_kd <- oe_vs_kd[,c('name', 'diff.log2.fc', 'diff.lg.p', 'diff.lg.fdr')]

met <- merge(kd_vs_ctrl, oe_vs_kd, by = 'name')
d <- duplicated(met$name)
met <- met[!d,]

table(met$diff.log2.fc.x*met$diff.log2.fc.y > 0)['TRUE']
table(met$diff.log2.fc.x*met$diff.log2.fc.y < 0)['TRUE']

met[met$diff.log2.fc.x*met$diff.log2.fc.y < 0,]

library(ensembldb)
library(EnsDb.Hsapiens.v86) 
# extract transcripts annotations
ensdb.hs <- EnsDb.Hsapiens.v86

# extract genes annotations
gene.hs <- genes(ensdb.hs, return.type = 'DataFrame')
gene.hs <- gene.hs[,c('gene_id', 'symbol', 'entrezid', 'gene_biotype')]

met <- merge(met, gene.hs, by.x = 'name', by.y = 'gene_id')
write.csv(met, 'm6a_compare_all.csv')

# clean mem
rm(gene.hs, ensdb.hs)
```
