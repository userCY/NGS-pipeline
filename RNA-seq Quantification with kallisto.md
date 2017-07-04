RNA-seq Quantification with kallisto
========================================

### 1. Building Kallisto human transcriptome index
```linux
kallisto index -i GRCh38_trans.idx <ref transcriptome>
```
### 2. Quantification with Kallisto
```R
kallisto quant -i <transcriptome index> -o output -b 100 --single -l 200 -s 20 -t 5 <input file>
```
Options:

  options  | meaning  
:---:|:-----:
-i <br> --index=STRING | Filename for the kallisto index to be used for pseudoalignment
-o <br> --output-dir=STRING | Directory to write output to
--single | Quantify single-end reads
-l <br> --fragment-length=DOUBLE | Estimated average fragment length
-s <br> --sd=DOUBLE | Estimated standard deviation of fragment length (default: -l, -s values are estimated <br> from paired end data, but are required when using --single)
-t <br> --threads=INT | Number of threads to use (default: 1)
