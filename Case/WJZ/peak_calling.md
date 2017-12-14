```R
library("exomePeak")

gtf=file.path('','home','ycli','data1','hisat2-2.1.0','indexes','Hsapiens_GRCh38_tran','Homo_sapiens.GRCh38.84.gtf')

ip1=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-01_ip.bam')
ip2=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-02_ip.bam')
ip3=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-03_ip.bam')

ip7=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-07_ip.bam')
ip8=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-08_ip.bam')
ip9=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-09_ip.bam')

in1=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-1_input.bam')
in2=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-2_input.bam')
in3=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-3_input.bam')

in4=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-4_input.bam')
in5=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-5_input.bam')
in6=file.path('','home','ycli','data1','WJZ','m6a_raw','clean','hisat2_out','analyze','wjz-6_input.bam')

result <- exomepeak(GENE_ANNO_GTF = gtf,
					IP_BAM = c(ip1,ip2,ip3),
                    INPUT_BAM = c(in1,in2,in3),
                    TREATED_IP_BAM = c(ip7,ip8,ip9),
                    TREATED_INPUT_BAM = c(in4,in5,in6),
                    OUTPUT_DIR = file.path('','home','ycli','data1','WJZ','kd_vs_ctrl_exomePeak_out'),
					SLIDING_STEP = 20,
					FRAGMENT_LENGTH = 150,
					READ_LENGTH = 150
                    )
save(result, file = file.path('','home','ycli','data1','WJZ','kd_vs_ctrl_exomePeak_out','result.rda'))

rmt <- RMT(IP_BAM = c(ip1,ip2,ip3,ip7,ip8,ip9),
			INPUT_BAM = c(in1,in2,in3,in4,in5,in6),
			GENE_ANNO_GTF = gtf,
			EXOME_OUTPUT_DIR = file.path('','home','ycli','data1','WJZ','exome_out')
			)
```


