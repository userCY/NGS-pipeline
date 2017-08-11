The differential expression analysis of an unconventional m6As-seq data
=================

### 1. read in kallisto tx expression data with tximport

```R
files <- file.path('D:','RWD', "yinrong", 'kallisto_output', 'count_matrix', row.names(sampleTable), "abundance_mod.tsv")
names(files) <- row.names(sampleTable)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene.mm)
head(txi.kallisto.tsv$counts)
```
### 2. gene annotation

using EnsDb.Mmusculus.v79 database and ensembldb package to annotate our genes
package manual: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

EnsDb gene features:

listTables(edb)

listColumns(edb, "gene")
```
$gene
[1] "gene_id"          "gene_name"        "entrezid"        
[4] "gene_biotype"     "gene_seq_start"   "gene_seq_end"    
[7] "seq_name"         "seq_strand"       "seq_coord_system"
[10] "symbol"
```

```R
library(EnsDb.Mmusculus.v79)

ensdb.mm <- EnsDb.Mmusculus.v79
tx.mm <- transcripts(ensdb.mm,return.type = "DataFrame")
tx.mm <- tx.mm[,c("tx_id", "gene_id")]

# extract genes annotations
gene.mm <- genes(ensdb.mm, return.type = 'DataFrame')
gene.mm <- gene.mm[,c('gene_id', 'symbol', 'entrezid', 'gene_biotype')]

# merge tx and gene annotations by gene_id
tx2gene.mm <- merge(tx.mm, gene.mm, by = 'gene_id')
tx2gene.mm <- tx2gene.mm[,c('tx_id', 'gene_id', 'symbol', 'entrezid', 'gene_biotype')]

# clean mem
rm(tx.mm, gene.mm, tx2gene.mm, ensdb.mm)

d <- duplicated(tx2gene.mm$gene_id)
tx2gene.mm <- tx2gene.mm[!d,]
```
### 3. extracting raw counts

raw\_count: raw read counts at gene level

0 filtered count: filter out 0 expressing genes across all samples

annotations: 
- ensembl gene id
- gene symbol
- Entrez gene id
- gene biological type

##### 3.1 input data
```R
raw_count <- txi.kallisto.tsv$counts
raw_count <- merge(raw_count, tx2gene.mm, by.x = 0, by.y = 'gene_id')
rn <- c(' ',' ',' ',' ','GMP-1', 'LSK-1', 'MEP-1', 'CMP-1', 'GMP-2', 'LSK-2', 'MEP-2', 'CMP-2', 'LT-HSC-1', 'PROG-1',
        'LIN-1','LT-HSC-2', 'PROG-2', 'LIN-2', 'LT-HSC-3', 'PROG-3', 'LIN-3', 'GMP-3', 'LSK-3', 'MEP-3',
        'CMP-3')

raw_count.input <- cbind(raw_count[,1], raw_count[,45:47],raw_count[2:22])
raw_count.input <- rbind(rn, raw_count.input)
raw_count.input[,5:25] <- raw_count.input[,5:25][,c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.input)[5:25] <- colnames(raw_count.input)[5:25][c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.input)[1] <- 'gene_id'
write.csv(raw_count.input, 'raw_count_input.csv')
rm(raw_count.input)

raw_count.input.clean <- cbind(raw_count[,1], raw_count[,45:47],raw_count[2:22])
raw_count.input.clean[,5:25] <- raw_count.input.clean[,5:25][,c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.input.clean)[5:25] <- colnames(raw_count.input.clean)[5:25][c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
raw_count.input.clean <- raw_count.input.clean[rowSums(raw_count.input.clean[5:25])>0,]
write.csv(raw_count.input.clean, 'raw_count_input_clean.csv')
rm(raw_count.input.clean)
```
##### 3.2 ip data

```R
raw_count.ip <- cbind(raw_count[,1], raw_count[,45:47],raw_count[23:43])
raw_count.ip <- rbind(rn, raw_count.ip)
raw_count.ip[,5:25] <- raw_count.ip[,5:25][,c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.ip)[5:25] <- colnames(raw_count.ip)[5:25][c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.ip)[1] <- 'gene_id'
write.csv(raw_count.ip, 'raw_count_ip.csv')
rm(raw_count.ip)

raw_count.ip.clean <- cbind(raw_count[,1], raw_count[,45:47],raw_count[23:43])
raw_count.ip.clean[,5:25] <- raw_count.ip.clean[,5:25][,c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
colnames(raw_count.ip.clean)[5:25] <- colnames(raw_count.ip.clean)[5:25][c(9,12,15,2,6,19,10,13,16,4,8,21,1,5,18,3,7,20,11,14,17)]
raw_count.ip.clean <- raw_count.ip.clean[rowSums(raw_count.ip.clean[5:25])>0,]
write.csv(raw_count.ip.clean, 'raw_count_ip_clean.csv')
rm(raw_count.ip.clean)
rm(raw_count)

raw_tpm <- txi.kallisto.tsv$abundance
rm(raw_tpm)
```

### 4. generating DGEList object for EdgeR

```R
library(edgeR)

# generate count matrix with gene annotation
cts.input <- txi.kallisto.tsv$counts[,1:21]
cts.input <- merge(cts.input, tx2gene.mm, by.x = 0, by.y = 'gene_id', all.x = TRUE, all.y = FALSE)
row.names(cts.input) <- cts.input[,1]
colnames(cts.input)
gene_anno.input <- cts.input[,c(1,24,25,26)]
colnames(gene_anno.input)[1] <- 'gene_id'
cts.input <- cts.input[,2:22]
colnames(cts.input) <- c('GMP-1', 'LSK-1', 'MEP-1', 'CMP-1', 'GMP-2', 'LSK-2', 'MEP-2', 'CMP-2', 'LT-HSC-1', 'PROG-1',
                         'LIN-1','LT-HSC-2', 'PROG-2', 'LIN-2', 'LT-HSC-3', 'PROG-3', 'LIN-3', 'GMP-3', 'LSK-3', 'MEP-3',
                         'CMP-3')

# normalize for gene length, seq depth, and RNA compensation
normMat <- txi.kallisto.tsv$length[,1:21]
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts.input/normMat)) + log(colSums(cts.input/normMat))

# biological gropu info
group.input <- c(as.character(sampleTable$condition[1:21]))

# create DGEList object
y.input <- DGEList(counts=cts.input,group=group.input, genes=gene_anno.input)

# feed calculated offset
y.input$offset <- t(t(log(normMat)) + o)

# filter out low-expressing genes (according to biological replicates)
# counts per million reads > 1 in more than 3 samples
keep <- rowSums(cpm(y.input)>1) > 3
y.input <- y.input[keep, , keep.lib.sizes=FALSE]
normMat <- normMat[keep,]
rm(keep)

# recompute the offsets (optional):
o1 <- log(calcNormFactors(y.input$counts/normMat)) + log(colSums(y.input$counts/normMat))
y.input$offset <- t(t(log(normMat)) + o1)

#y.input.cpm <- cpm(y.input, normalized.lib.sizes = T)
#y.input.rpkm <- rpkm(y.input, normalized.lib.sizes = T, gene.length = y.input$genes)
#write.csv(y.input.rpkm, 'input_rpkm.csv')

# plot PCA
plotMDS(y.input)
```

