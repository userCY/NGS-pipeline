import os

os.system('mkdir deep_out')

groups = ['zt-07.bam zt-14.bam zt-16.bam', 'zt-18.bam zt-02.bam zt-04.bam', 'zt-06.bam zt-10.bam zt-12.bam', 
'zt-08.bam zt-15.bam zt-17.bam', 'zt-01.bam zt-03.bam zt-05.bam', 'zt-09.bam zt-11.bam zt-13.bam']

labels = ['C1-input C1-F-IP C1-M-IP', 'C2-input C2-F-IP C2-M-IP', 'C3-input C3-F-IP C3-M-IP' , 
'P1-input P1-F-IP P1-M-IP', 'P2-input P2-F-IP P2-M-IP', 'P3-input P3-F-IP P3-M-IP']
plot_names = ['C1','C2','C3','P1','P2','P3']
n = len(groups)

for i in range(n):
	cmd = '/home/ycli/data1/deepTools/bin/plotFingerprint \
	-p 20 \
	-b {0} \
	--labels {1} \
	--minMappingQuality 30 --skipZeros \
	--region 19 --numberOfSamples 50000 \
	-T "Fingerprints of different samples"  \
	--plotFile ./deep_out/fingerprints-{2}.png'.format(groups[i], labels[i], plot_names[i])
	print('ploting {0} for {1}'.format(groups[i], labels[i]))
	os.system(cmd)
