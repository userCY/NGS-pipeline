kallisto2edgeR with tximport
====================

### 1. import required libraries

```R
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86) 
# choose the appropriate or latest ensembl annotation
```
### 2. create tx to gene symbol annotations from ensembl database

```R
# extract transcripts annotations
ensdb.hs <- EnsDb.Hsapiens.v86
tx.hs <- transcripts(ensdb.hs,return.type = "DataFrame")
tx.hs <- k[,c("tx_id", "gene_id")]

# extract genes annotations
gene.hs <- genes(ensdb.hs, return.type = 'DataFrame')
gene.hs <- gene.hs[,c('gene_id', 'symbol')]

# merge tx and gene annotations by gene_id
tx2gene.hs <- merge(tx.hs, gene.hs, by = 'gene_id')

# only keep tx ID, gene ID
tx2gene <- tx2genes.hs[, 2:3]  # tx ID, then gene ID

# clean mem
rm(tx.hs, gene.hs, tx2gene.hs)
```

### :exclamation: 3. debug: solve the problem of ENST + version bumber, remove version number
> **Explanation:**

> Somewhere in the last two years, Ensembl and Gencode (and I think NCBI is doing this now too) decided to release FASTA and GTF files that include both the Ensembl ID and the version number (eg. "ENSMUST00000000001.4" instead of just "ENSMUST00000000001"). 

> The IDs with both the ID number and the version number are what's used in Kallisto, and will be your "target_ids".

The kallisto quantification tx_id doesn't match the tx2gene map tx_id, which will casue an error. 
The problem is to be addresed with:
```R
abund_file <- read.table(file = files[4], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
abund_file[,1] <- substr(abund_file[,1], 1, 15)
write.table(abund_file, file='abundance.tsv', quote=FALSE, sep='\t', row.names = FALSE)
```
**_Note that these codes create an abundance.tsv file under the R working directory.
The modified file need to be placed back to the kallisto output directory manually_**

### 4. creating txi file from kallisto output directory
```R
files <- file.path('D:','RWD', "FTO", row.names(sampleTable), "abundance.tsv")
names(files) <- row.names(sampleTable)
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene)
head(txi.kallisto.tsv$counts)
```

sampleTable file format:
--------------------------
SampleTable data.frame is used in kallisto, sleuth and tximport.
The basic format of this df contains only 1 column:
rownames: sample name
colnames: 'condition'
colData: condition group relating to the sample (eg. control, KO)

FTO example:
|  | condition |
| ------------- | ------------- |
| FTO+\_1 | FTO+ |
| FTO+\_2 | FTO+ |
| WT_1 | WT |
| WT_2 | WT |
