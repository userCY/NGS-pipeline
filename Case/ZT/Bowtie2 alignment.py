```python
import os 
import re
import fnmatch

files = os.listdir(os.getcwd())
forward = fnmatch.filter(files, '*_R1.fastq.gz.paired')
forward.sort()
reverse = fnmatch.filter(files, '*_R2.fastq.gz.paired')
reverse.sort()

bowtie2_idx = '/home/ycli/data1/bowtie2-2.3.2/index/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
output_dir_base = 'bowtie2_out'
n = len(forward)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{2}-\d{2}')
	m = p.search(forward[i])
	outfile = m.group().lower()
	bowtie2_cmd = 'bowtie2 ' + '-p 20 -x ' + bowtie2_idx + ' -1 ' + forward[i] + ' -2 ' + reverse[i] + ' -S ./' + output_dir_base + '/' + outfile + '.sam' + ' 2>' + './' + output_dir_base + '/' + outfile + 'align.log'
	os.system(bowtie2_cmd)
  
  ```
