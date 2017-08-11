The differential expression analysis of an unconventional m6As-seq data
=================

1. read in kallisto tx expression data with tximport

```R
files <- file.path('D:','RWD', "yinrong", 'kallisto_output', 'count_matrix', row.names(sampleTable), "abundance_mod.tsv")
names(files) <- row.names(sampleTable)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene.mm)
head(txi.kallisto.tsv$counts)
```
2. gene annotation

using EnsDb.Mmusculus.v79 data base and ensembldb package to annotate our genes
package manual: https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

EnsDb gene features:

$gene

[1] "gene_id"          "gene_name"        "entrezid"        
[4] "gene_biotype"     "gene_seq_start"   "gene_seq_end"    
[7] "seq_name"         "seq_strand"       "seq_coord_system"
[10] "symbol"

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
