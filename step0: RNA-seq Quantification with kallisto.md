RNA-seq Quantification with kallisto
========================================

### 1. Building Kallisto human transcriptome index
```R
kallisto index -i GRCh38_trans.idx <ref transcriptome>
```
### 2. Quantification with Kallisto
for single-end sequencing:
```R
kallisto quant -i <transcriptome index> -o output -b 100 --single -l 200 -s 20 -t 5 <input file>
```
for paired-end sequencing:
```R
kallisto quant -i <transcriptome index> -o output -b 10 -t 5 <paired_1> <paired_2>
```
**note: -l -s options are not needed in paired-end mode 'cause the programme will automatically estimates their values**

Options:

  options  | meaning  
:---:|:-----:
-i <br> --index=STRING | Filename for the kallisto index to be used for pseudoalignment
-o <br> --output-dir=STRING | Directory to write output to
--single | Quantify single-end reads
-l <br> --fragment-length=DOUBLE | Estimated average fragment length
-s <br> --sd=DOUBLE | Estimated standard deviation of fragment length (default: -l, -s values are estimated <br> from paired end data, but are required when using --single)
-t <br> --threads=INT | Number of threads to use (default: 1)

### 3. For batch processing
***Example:***
```python
# wording directory: fastq file directory

import os
import re
import fnmatch

# extract forward and reverse seq file names
file = os.listdir(os.getcwd())
foward = fnmatch.filter(file, '*_1.fq.gz')
reverse = fnmatch.filter(file, '*_2.fq.gz')

# determine file number
n = len(forward)

# set file path
kallisto_idx = 'index path'
output_dir_base = 'out dir base name'


# batch processing
for i in range(n)：
	p = re.compile(r'\w\w\d{1,2}')
	m = p.search(forward[i])
	outdir = m.group().lower()
	kallisto_cmd = 'kallisto quant -i'+' '+'' -o'+' '+kallisto_idx + output_dir_base + outdir+' '+'-b 10 -t 10'
	os.system(kallisto_cmd+' '+forward[i]+' '+reverse[i])
```

kallisto output file format:
----------------
kallisto quant produces three output files by default:

- **abundances.h5** is a HDF5 binary file containing run info, abundance esimates, bootstrap estimates, and transcript length information length. This file can be read in by sleuth
- **abundances.tsv** is a plaintext file of the abundance estimates. It does not contains bootstrap estimates. Please use the --plaintext mode to output plaintext abundance estimates. Alternatively, kallisto h5dump can be used to output an HDF5 file to plaintext. The first line contains a header for each column, including estimated counts, TPM, effective length.
- **run_info.json** is a json file containing information about the run

**note: either abundances.h5 or abundances.tsv can be used to construct count matrix in subsequent analysis**
