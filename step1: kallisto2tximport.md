kallisto output data processed by tximport for downstream analysis
====================

### 1. import required libraries

```R
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86) 
```
>The tximport package has a single function for importing transcript-level estimates. The type argument is used to specify what software was used for estimation (“kallisto”, “salmon”, “sailfish”, and “rsem” are implemented). A simple list with matrices, “abundance”, “counts”, and “length”, is returned, where the transcript level information is summarized to the gene-level. The “length” matrix can be used to generate an offset matrix for downstream gene-level differential analysis of count matrices, as shown below.

**note:**
> tximport could alternatively generate counts from abundances, using the argument countsFromAbundance, scaled to library size, "scaledTPM", or additionally scaled using the average transcript length, averaged over samples and to library size, "lengthScaledTPM". Using either of these approaches, **the counts are not correlated with length**, and so the length matrix should not be provided as an offset for downstream analysis packages. 

### 2. create tx to gene symbol annotations from ensembl database

```R
# extract transcripts annotations
ensdb.hs <- EnsDb.Hsapiens.v86
tx.hs <- transcripts(ensdb.hs,return.type = "DataFrame")
tx.hs <- tx.hs[,c("tx_id", "gene_id")]

# extract genes annotations
gene.hs <- genes(ensdb.hs, return.type = 'DataFrame')
gene.hs <- gene.hs[,c('gene_id', 'symbol')]

# merge tx and gene annotations by gene_id
tx2gene.hs <- merge(tx.hs, gene.hs, by = 'gene_id')

# only keep tx ID, gene ID
tx2gene <- tx2gene.hs[, 2:3]  # tx ID, then gene ID

# clean mem
rm(tx.hs, gene.hs, tx2gene.hs)
```

### :exclamation: 3. debug: solve the problem of ENST + version bumber, remove version number
> **Explanation:**

> Somewhere in the last two years, Ensembl and Gencode (and I think NCBI is doing this now too) decided to release FASTA and GTF files that include both the Ensembl ID and the version number (eg. "ENSMUST00000000001.4" instead of just "ENSMUST00000000001"). 

> The IDs with both the ID number and the version number are what's used in Kallisto, and will be your "target_ids".

The kallisto quantification tx_id doesn't match the tx2gene map tx_id, which will casue an error. 
The problem is to be addresed with:

create a sample table first (format: rownames:sample; colnames:condition; colData: group) and corresponding file path:
```R
sampleTable <- data.frame('condition' = c('FTO+', 'FTO+', 'WT', 'WT'), stringsAsFactors = FALSE)
row.names(sampleTable) <- c('FTO+_1', 'FTO+_2', 'WT_1', 'WT_2')
file <- file.path('D:', 'RWD', row.names(sampleTable), 'abundance.tsv')
```
removing version number:
```R
abund_file <- read.table(file = file[4], sep = '\t', header = TRUE, stringsAsFactors = FALSE)
abund_file[,1] <- substr(abund_file[,1], 1, 15)
write.table(abund_file, file='abundance.tsv', quote=FALSE, sep='\t', row.names = FALSE)
```
**_Note that these codes create an abundance.tsv file under the R working directory.
The modified file need to be placed back to the kallisto output directory manually_**

### 4. creating txi file from kallisto output directory
```R
files <- file.path('D:','RWD', "FTO", row.names(sampleTable), "abundance.tsv")
names(files) <- row.names(sampleTable)
```
generate raw counts with:
```R
txi.kallisto.tsv <- tximport(file, type = "kallisto", tx2gene = tx2gene)
head(txi.kallisto.tsv$counts)
```
generate length-independent scaled counts with:
```R
txi <- tximport(file, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
head(txi.kallisto.tsv$counts)
```
scaled counts are suitable for gene expression related plots, such as heatmap or expression bar plot

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
