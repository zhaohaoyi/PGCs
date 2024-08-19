library(ggplot2)
library(ggpubr)

data<-read.table("exp.xls", sep="\t", header=T, check=F, comment="", na.strings="", stringsAsFactors=F)

p = ggviolin(data, x = "", y = "Expression Level(log10)", 
			 fill = "Cluster", color = "Cluster", alpha = 1, width = 0.75) 
ggsave(p, filename = "violin.pdf", height = 6, width=8)

