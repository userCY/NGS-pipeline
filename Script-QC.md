Python Script - FastQC
==================

```python
import os

# in fa file folder
# get fa file list
files = os.listdir(os.getcwd())

# make QC output file
os.system('mkdir QC_output')

for f in files:
  os.system('fastqc '+ f + ' -o QC_output')
```

Python Script - Trimmomatic
===========================

```python
import os 
import re

# adapter
adapter = '/home/ycli/data1/Trimmomatic-0.36/adapters/nextra-primer.fa'

#headcrop
headcrop = 15

#fq files
files = os.listdir(os.getcwd())
fqfiles = []
p = re.compile('\w+\.fastq\.gz')
for f in files:
	result = p.search(f)
	if result != None:
		fqfiles.append(result.group())
	
#forward fq files
p = re.compile('\w+_R1.\w+.\w+')
forward = []
for f in fqfiles:
	result = p.search(f)
	if result != None:
		forward.append(result.group())
	
# reverse fq files
p = re.compile('\w+_R2.\w+.\w+')
reverse = []
for f in fqfiles:
	result = p.search(f)
	if result != None:
		reverse.append(result.group())

# sort
forward.sort()
reverse.sort()


# linux command
os.system('mkdir clean')
for i,f in enumerate(forward):
	fq1 = f
	print fq1
	print i
	fq2 = reverse[i]
	out1 = './clean/'+fq1+'.paired'
	out12 = './clean/'+fq1+'.unpaired'
	out2 = './clean/'+fq2+'.paired'
	out22 = './clean/'+fq2+'.unpaired'
	print 'processing trimming {0} {1}'.format(fq1,fq2)
	print
	cmd='java -jar /home/ycli/data1/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads 10\
	{0} {1} {2} {3} {4} {5} \
	HEADCROP:{6} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 ILLUMINACLIP:{7}:2:30:10'.format(fq1,fq2,out1,out12,out2,out22,headcrop,adapter)
	os.system(cmd)
```
