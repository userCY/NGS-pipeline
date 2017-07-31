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
