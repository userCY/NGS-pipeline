Python script RNA-seq with Kallisto
===============================


```python
import os 
import re
import fnmatch

files = os.listdir(os.getcwd())
forward = fnmatch.filter(files, '*_R1.fastq.gz.paired')
forward.sort()
reverse = fnmatch.filter(files, '*_R2.fastq.gz.paired')
reverse.sort()

kallisto_idx = '/home/ycli/data1/kallisto_linux-v0.43.1/GRCh38_trans.idx'
output_dir_base = 'kallisto_out'
n = len(forward)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{3}-\d{1,2}')
	m = p.search(forward[i])
	outdir = m.group().lower()
	kallisto_cmd = 'kallisto quant -i' + kallisto_idx + ' ' + '-o' + ' ./' + output_dir_base + '/' + outdir +' '+'-b 10 -t 10'
	os.system(kallisto_cmd+' '+forward[i]+' '+reverse[i])

```
