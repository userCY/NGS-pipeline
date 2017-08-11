The differential expression analysis of an unconventional m6As-seq data
=================

1. read in kallisto tx expression data with tximport

```python
files <- file.path('D:','RWD', "yinrong", 'kallisto_output', 'count_matrix', row.names(sampleTable), "abundance_mod.tsv")
names(files) <- row.names(sampleTable)

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene.mm)
head(txi.kallisto.tsv$counts)
```
2. gene annotation

```python


d <- duplicated(tx2gene.mm$gene_id)
tx2gene.mm <- tx2gene.mm[!d,]

```
