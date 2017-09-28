
import os 
import re

#fq files
files = os.listdir(os.getcwd())
p1 = re.compile('\w+-\d+.sam')
p2 = re.compile('\w+-\d+.')
samfiles = []
for f in files:
	result1 = p1.search(f)
	if result1 != None:
		samfiles.append(result1.group())
outfiles = []
for f in files:
	result2 = p2.search(f)
	if result2 != None:
		outfiles.append(result2.group())
		
# sort
samfiles.sort()
outfiles.sort()

# linux command
file_path = os.getcwd()
for i,f in enumerate(samfiles):
	print 'compressing {0}'.format(f)
	print
	cmd = 'samtools sort -o {0}/{2}bam -@ 20 {0}/{1}'.format(file_path,f,outfiles[i])
	os.system(cmd)


#fq files
files = os.listdir(os.getcwd())
p3 = re.compile('\w+-\d+.bam')
#p2 = re.compile('\w+-\d+_\w+.')
bamfiles = []
for f in files:
	result1 = p3.search(f)
	if result1 != None:
		bamfiles.append(result1.group())

# linux command
file_path = os.getcwd()
for f in bamfiles:
	print 'indexing {0}'.format(f)
	print
	cmd = 'samtools index {0}/{1} -@ 20'.format(file_path,f)
	os.system(cmd)
