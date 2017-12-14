This script plot sequencing depth form bedgraph files using R package ggplot2
======

```R
library(ggplot2)

# chr start end depth

rgb(red=230, green=75, blue=52, maxColorValue = 255)
rgb(red=60, green=83, blue=136, maxColorValue = 255)

ip_kd <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-09_ip.bedgraph', sep = '\t')
ip_kd_19 <- ip_kd[ip_kd$V1 == 19,]
tail(ip_kd_19[ip_kd_19$V2 <= 41219203,])
head(ip_kd_19[ip_kd_19$V3 >= 41261766,])
ip_kd_axl <- ip_kd_19[ip_kd_19$V2 > 41217870 & ip_kd_19$V3 < 41262860,]

in_kd <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-6_input.bedgraph', sep = '\t')
in_kd_19 <- in_kd[in_kd$V1 == 19,]
tail(in_kd_19[in_kd_19$V2 <= 41219203,])
head(in_kd_19[in_kd_19$V3 >= 41261766,])
in_kd_axl <- in_kd_19[in_kd_19$V2 > 41211195 & in_kd_19$V3 < 41262585,]

ggplot() +
  geom_step(data = ip_kd_axl, mapping = aes(x = V2, y = V4), color = '#E64B34', size = 2) +
  geom_step(data = in_kd_axl, mapping = aes(x = V2, y = V4), color = '#3C5388', size = 2) +
  theme_classic() +
  xlim(41219203, 41261766) +
  ylim(0, 400)

ip_con <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-02_ip.bedgraph', sep = '\t')
ip_con_19 <- ip_con[ip_con$V1 == 19,]
tail(ip_con_19[ip_con_19$V2 <= 41219203,])
head(ip_con_19[ip_con_19$V3 >= 41261766,])
ip_con_axl <- ip_con_19[ip_con_19$V2 > 41207530 & ip_con_19$V3 < 41264325,]

in_con <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-3_input.bedgraph', sep = '\t')
in_con_19 <- in_con[in_con$V1 == 19,]
tail(in_con_19[in_con_19$V2 <= 41219203,])
head(in_con_19[in_con_19$V3 >= 41261766,])
in_con_axl <- in_con_19[in_con_19$V2 > 41215705 & in_con_19$V3 < 41262570,]

ggplot() +
  geom_step(data = ip_con_axl, mapping = aes(x = V2, y = V4), color = '#E64B34', size = 2) +
  geom_step(data = in_con_axl, mapping = aes(x = V2, y = V4), color = '#3C5388', size = 2) +
  theme_classic() +
  xlim(41219203, 41261766) +
  ylim(0, 400)

ip_res <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-05_ip.bedgraph', sep = '\t')
ip_res_19 <- ip_res[ip_res$V1 == 19,]
tail(ip_res_19[ip_res_19$V2 <= 41219203,])
head(ip_res_19[ip_res_19$V3 >= 41261766,])
ip_res_axl <- ip_res_19[ip_res_19$V2 > 41207525 & ip_res_19$V3 < 41264275,]

in_res <- read.table('C:\\Users\\wanvdphelys\\Desktop\\Zhanglab\\Wang Jiazhen\\bdg_out\\wjz-8_input.bedgraph', sep = '\t')
in_res_19 <- in_res[in_res$V1 == 19,]
tail(in_res_19[in_res_19$V2 <= 41219203,])
head(in_res_19[in_res_19$V3 >= 41261766,])
in_res_axl <- in_res_19[in_res_19$V2 > 41207510 & in_res_19$V3 < 41262585,]

ggplot() +
  geom_step(data = ip_res_axl, mapping = aes(x = V2, y = V4), color = '#E64B34', size = 2) +
  geom_step(data = in_res_axl, mapping = aes(x = V2, y = V4), color = '#3C5388', size = 2) +
  theme_classic() +
  xlim(41219203, 41261766) +
  ylim(0, 400)
```
