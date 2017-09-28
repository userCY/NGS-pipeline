import os 
import re
import fnmatch

files = os.listdir(os.getcwd())
bamfiles = fnmatch.filter(files, '*.bam')
bamfiles.sort()

output_dir_base = 'bigwig_out'
n = len(bamfiles)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{2}-\d{2}')
	m = p.search(bamfiles[i])
	outfile = m.group().lower()
	bamCoverage_cmd = 'bamCoverage -p 20 -b ' + bamfiles[i] + ' -o ' + './' + output_dir_base + '/' + outfile + '.bigwig' + ' --binSize 5 --normalizeUsingRPKM --outFileFormat bigwig'
	os.system(bamCoverage_cmd)
	
	
output_dir_base = 'bdg_out'
n = len(bamfiles)

os.system('mkdir ' + output_dir_base)
for i in range(n):
	p = re.compile(r'\w{2}-\d{2}')
	m = p.search(bamfiles[i])
	outfile = m.group().lower()
	bamCoverage_cmd = 'bamCoverage -p 20 -b ' + bamfiles[i] + ' -o ' + './' + output_dir_base + '/' + outfile + '.bedgraph' + ' --binSize 5 --normalizeUsingRPKM --outFileFormat bedgraph'
	os.system(bamCoverage_cmd)
