# -*- coding: utf-8 -*-
import os 
import re
import fnmatch

files = os.listdir(os.getcwd())
bwfiles = fnmatch.filter(files, '*.bigwig')
bwfiles.sort()

output_dir_base = 'deep_matrix'
n = len(bwfiles)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{2}-\d{2}')
	m = p.search(bwfiles[i])
	outfile = m.group().lower()
	computeMatrix_cmd = 'computeMatrix scale-regions \
		-p 30 \
		-R /home/ycli/data1/hisat2-2.1.0/indexes/Hsapiens_GRCh38_tran/Homo_sapiens.GRCh38.84.gtf \
		-S {0} \
		-b 3000 -a 3000 \
		--regionBodyLength 5000 \
		--binSize 5 \
		--skipZeros \
		-o ./{1}/{2}_scaled.gz \
		--outFileNameMatrix ./{1}/{2}_scaled.tab \
		--outFileSortedRegions ./{1}/regions2_{2}.bed'.format(bwfiles[i],output_dir_base,outfile)
	print('computing {0} with a bin size of 5'.format(bwfiles[i]))
	os.system(computeMatrix_cmd)
