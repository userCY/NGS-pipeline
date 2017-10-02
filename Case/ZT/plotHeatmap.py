# -*- coding: utf-8 -*-
import os 
import re
import fnmatch

files = os.listdir(os.getcwd())
mtfiles = fnmatch.filter(files, '*.gz')
mtfiles.sort()

output_dir_base = 'heatmap'
n = len(mtfiles)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{2}-\d{2}')
	m = p.search(mtfiles[i])
	outfile = m.group().lower()
	plotHeatmap_cmd = 'plotHeatmap -m {0} \
		-out ./{1}/{2}_Heatmap2.png \
		--colorMap RdBu \
		--zMin -3 --zMax 3 \
		--dpi 300'.format(mtfiles[i], output_dir_base, outfile)
	os.system(plotHeatmap_cmd)
