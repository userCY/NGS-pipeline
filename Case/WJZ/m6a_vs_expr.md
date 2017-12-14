Plot mRNA m6a level against gene expression level

1 read in kallisto file
==================

```R
sampleTable <- data.frame('condition' = c('control', 'control', 'KD', 'KD'), stringsAsFactors = FALSE)
row.names(sampleTable) <- c('ctrl_1', 'ctrl_2', 'KD_1', 'KD_2')
file <- file.path('D:', 'RWD', 'WJZ','wjz_kallisto',c('control-2.tsv', 'control-3.tsv','KD-1.tsv','KD-3.tsv'))
names(file) <- row.names(sampleTable)

d <- duplicated(tx2gene.hs$gene_id)
tx2gene.hs <- tx2gene.hs[!d,]

txi.kallisto.tsv <- tximport(file, type = "kallisto", tx2gene = tx2gene.hs)
```

2 differential expression analysis with edgeR
============

```R
raw_count <- txi.kallisto.tsv$counts

raw_count <- merge(raw_count, tx2gene.hs, by.x = 0, by.y = 'gene_id')

cts <- raw_count[2:5]
row.names(cts) <- raw_count[,1]
cts <- merge(cts, tx2gene.hs, by.x = 0, by.y = 'gene_id', all.x = TRUE, all.y = FALSE)
gene_anno <- cts[,c(1,7,8,9)]
colnames(gene_anno)[1] <- 'gene_id'
cts <- cts[,2:5]

normMat <- txi.kallisto.tsv$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(as.matrix(cts)/normMat)) + log(colSums(as.matrix(cts)/normMat))

group <- c(as.character(sampleTable$condition))

y <- DGEList(counts=cts,group=group, genes=gene_anno)

y$offset <- t(t(log(normMat)) + o)

keep <- rowSums(cpm(y)>1) > 2
y <- y[keep, , keep.lib.sizes=FALSE]
normMat <- normMat[keep,]

o1 <- log(calcNormFactors(y$counts/normMat)) + log(colSums(y$counts/normMat))
y$offset <- t(t(log(normMat)) + o1)

plotMDS(y)

design <- model.matrix(~0+group)

y <- estimateDisp(y, design)
y$common.dispersion

sqrt(y$common.disp)
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast = c(-1,1))

plotMD(lrt, main = 'knockdown vs control')
abline(h=c(-1, 1), col="blue")

# count upregulated (1) and downregulated (-1) gene numbers in KD vs ctrl
summary(decideTests(lrt))
```
3 starburst plot
============

```R
exp <- topTags(lrt, 10000000)$table
met <- read.table('D:\\RWD\\WJZ\\diff_peak.xls', header = T, stringsAsFactors = F, sep = '\t')
met <- met[,c('name', 'diff.log2.fc', 'diff.lg.p', 'diff.lg.fdr')]
met$diff.lg.p <- 10^met$diff.lg.p
met$diff.lg.fdr <- 10^met$diff.lg.fdr
colnames(met) <- c('gene_id', 'logFC', 'PValue', 'FDR')


# define parameters

.e <- environment()

group1 <- 'control'
group2 <- 'KD'


exp.p.cut <- 0.05
met.p.cut <- 0.05

# diffmean.cut <- 0
logFC.me.cut <- 1.5
#logFC.cut
logFC.ge.cut <- 1.5
# add names
names <- FALSE
names.fill = TRUE
# no circle
circle = FALSE
filename <- "starburst.png"
return.plot <- FALSE
ylab <- expression(atop("Gene Expression",
                       paste(Log[10],
                             " (FDR corrected P values)")))
xlab = expression(atop("mRNA adenine N6-methylation level",
                       paste(Log[10],
                             " (FDR corrected P values)")))
height<-10
width<-15
dpi<-300

title <- "ALKBH5 knockdown vs control Starburst Plot"
legend <- "m6a level/Expression Relation"

label <- c("Not Significant",
          "Up regulated & Hypo methylated",
          "Down regulated & Hypo methylated",
          "hypo methylated",
          "hyper methylated",
          "Up regulated",
          "Down regulated",
          "Up regulated & Hyper methylated",
          "Down regulated & Hyper methylated")

color <- c("#000000", "#E69F00","#56B4E9", "#009E73",
            "red", "#0072B2","#D55E00", "#CC79A7",
            "purple")

names(color) <- as.character(1:9)
names(label) <- as.character(1:9)
names.color <- color
names(names.color) <- label
#---------------------------------------------

volcano <- merge(met, exp, by = 'gene_id')
volcano <- volcano[!duplicated(volcano$symbol),]
#volcano$ID <- paste(volcano$Gene_symbol,
#                    volcano$probeID, sep = ".")


volcano$geFDR <- log10(volcano$FDR.y)
volcano$geFDR2 <- volcano$geFDR
volcano[volcano$logFC.y > 0, "geFDR2"] <- -1 * volcano[volcano$logFC.y > 0, "geFDR"]

volcano$meFDR <- log10(volcano$FDR.x)
volcano$meFDR2 <- volcano$meFDR
volcano[volcano$logFC.x > 0, "meFDR2"] <- -1 * volcano[volcano$logFC.x > 0, "meFDR"]

label[2:9] <-  paste(label[2:9], "in", group2)

met.lowerthr <- log10(met.p.cut)
met.upperthr <- (-met.lowerthr)

exp.lowerthr <- log10(exp.p.cut)
exp.upperthr <- (-exp.lowerthr)

# Group 2:up regulated and hypomethylated
a <- subset(volcano,
            volcano$geFDR2 > exp.upperthr &
              volcano$meFDR2 < met.lowerthr)

a.sig <- subset(a, abs(a$logFC.x) > logFC.me.cut &
                  abs(a$logFC.y) > logFC.ge.cut)


# Group 3: down regulated and hypomethylated
b <- subset(volcano,
            volcano$geFDR2 < exp.lowerthr &
              volcano$meFDR2 < met.lowerthr)

b.sig <- subset(b, abs(b$logFC.x) > logFC.me.cut &
                  abs(b$logFC.y) > logFC.ge.cut)


# Group 4: hypomethylated
c <- subset(volcano,
            volcano$geFDR2 > exp.lowerthr &
              volcano$geFDR2 < exp.upperthr &
              volcano$meFDR2 < met.lowerthr)

# Group 5: hypermethylated
d <- subset(volcano,
            volcano$geFDR2 > exp.lowerthr &
              volcano$geFDR2 < exp.upperthr &
              volcano$meFDR2 > met.upperthr)

# Group 6: upregulated
e <- subset(volcano,
            volcano$geFDR2 > exp.upperthr &
              volcano$meFDR2 < met.upperthr &
              volcano$meFDR2 > met.lowerthr)

# Group 7: downregulated
f <- subset(volcano,
            volcano$geFDR2 < exp.lowerthr &
              volcano$meFDR2 < met.upperthr &
              volcano$meFDR2 > met.lowerthr)

# Group 8: upregulated and hypermethylated
g <- subset(volcano,
            volcano$geFDR2 > exp.upperthr &
              volcano$meFDR2 > met.upperthr)

g.sig <- subset(g, abs(g$logFC.x) > logFC.me.cut &
                  abs(g$logFC.y) > logFC.ge.cut)

# Group 9: downregulated and hypermethylated
h <- subset(volcano,
            volcano$geFDR2 < exp.lowerthr &
              volcano$meFDR2 > met.upperthr)

h.sig <- subset(h, abs(h$logFC.x) > logFC.me.cut &
                  abs(h$logFC.y) > logFC.ge.cut)

groups <- as.character(seq(2,9))

# return methylation < 0, expressao >0
volcano$starburst.status  <-  "Not Significant"
volcano$shape <- "1"
volcano$threshold.starburst <- "1"
volcano$threshold.size <- "1"

state <- c("Hypo methylated & Up regulated",
           "Hypo methylated & Down regulated",
           "hypo methylated",
           "hyper methylated",
           "Up regulated",
           "Down regulated",
           "Hyper methylated & Up regulated",
           "Hyper methylated & Down regulated")

s <- list(a, b, c, d, e, f, g, h)
for (i in seq_along(s)) {
  idx <- rownames(s[[i]])
  if (length(idx) > 0) {
    volcano[idx, "threshold.starburst"] <- groups[i]
    volcano[idx, "starburst.status"] <-  state[i]
  }
}

size <- rep(2,4)
shape <-  as.character(rep(2,4))
volcano_normal <- volcano
significant <- NULL
s <- list(a.sig,b.sig,g.sig,h.sig)
for (i in seq_along(s)) {
  idx <- rownames(s[[i]])
  if (length(idx) > 0) {
    volcano[idx, "threshold.size"] <- size[i]
    volcano[idx, "shape"] <-  shape[i]
    significant <-  rbind(significant,volcano[idx,])
  }
}

library(ggplot2)
## starburst plot
p <- ggplot(data = volcano_normal, environment = .e,
            aes(x = volcano_normal$meFDR2,
                y = volcano_normal$geFDR2,
                colour = volcano_normal$threshold.starburst)) +
  geom_point()
#p <- p + scale_shape_discrete(
#    labels = c("Candidate Biologically Significant"),
#    name = "Biological importance")

if(!is.null(significant) & circle){
  p <- p + geom_point( data = significant,
                       aes(x = significant$meFDR2,
                           y = significant$geFDR2),
                       color = "black",
                       shape=1,
                       size = 8,
                       show.legend = FALSE)
}


library(ggrepel)
if(names == TRUE & !is.null(significant)){
  message("Adding names to genes")
  if(names.fill){
    p <- p + geom_label_repel(
      data = significant,
      aes(x = significant$meFDR2, y =  significant$geFDR2,
          label = significant$symbol, fill = as.factor(significant$starburst.status)),
      size = 4, show.legend = FALSE,
      fontface = 'bold', color = 'black',
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    ) + scale_fill_manual(values=names.color)
  }  else {
    p <- p + geom_text_repel(
      data = significant,
      aes(x = significant$meFDR2, y =  significant$geFDR2,
          label = significant$symbol, fill = significant$starburst.status),
      size = 4, show.legend = FALSE,
      fontface = 'bold', color = 'black',
      point.padding = unit(0.3, "lines"),
      box.padding = unit(0.5, 'lines')
    )
  }
}

if (!is.null(xlim)) {
  p <- p + xlim(xlim)
}
if (!is.null(ylim)) {
  p <- p + ylim(ylim)
}
p <- p + ggtitle(title) + ylab(ylab) + xlab(xlab) + guides(size=FALSE)
p <- p + scale_color_manual(values = color, labels = label, name = legend) +
  guides(col = guide_legend(nrow = 3))


p <-  p + geom_hline(aes(yintercept = exp.lowerthr), colour = "black",
                     linetype = "dashed") +
  geom_hline(aes(yintercept = exp.upperthr), colour = "black",
             linetype = "dashed") +
  geom_vline(aes(xintercept = met.lowerthr), colour = "black",
             linetype = "dashed") +
  geom_vline(aes(xintercept = met.upperthr), colour = "black",
             linetype = "dashed")  +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x=element_line(colour = "black"),
        axis.line.y=element_line(colour = "black"),
        legend.position="top",
        legend.key = element_rect(colour = 'white'),
        plot.title = element_text(face = "bold", size = 16,hjust = 0.5),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text= element_text(size = 14),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.text.x = element_text(vjust = 0.5,
                                   size = 14),
        axis.title.y = element_text(face = "bold",
                                    size = 14),
        axis.text.y = element_text(size = 14))

if(!return.plot) ggsave(filename = filename, width = width, height = height, dpi = dpi)

volcano$shape <- NULL
volcano$threshold.starburst <- NULL
volcano$threshold.size <- NULL

volcano <- subset(volcano, volcano$geFDR <= exp.lowerthr &
                    volcano$meFDR <= met.lowerthr)

if (logFC.me.cut != 0) {
  volcano <- subset(volcano, abs(volcano$logFC.x) > logFC.me.cut)
}
if (logFC.ge.cut != 0){
  volcano <- subset(volcano, abs(volcano$logFC.y) >= logFC.ge.cut)
}

if(return.plot) {
  return(list(plot=p,starburst=volcano))
}
```
