library(dplyr)
library(ggplot2)

file <- "data.xls"
outpre <- "UMAP"

data <- read.table(file,head =T,sep="\t",quote="",stringsAsFactors=T,check=F)
colnames(data)[2:3] <- c("x","y")
data$Phase <- factor(data$Phase, levels=c("G1","S","G2","M", "non-cycling"))

my_theme <- theme_bw() +
	theme( panel.background = element_blank(),
		legend.title =  element_blank() ,
		panel.grid = element_blank(),
		strip.background = element_blank(),
		strip.text =  element_text(size=13)
		)

colors = c("#00468b", "#42b540", "#ede447", "#ff7777", "#bebebe")
names(colors) = levels(data$Phase)
p <- ggplot(data) + 
	geom_point(aes(x = x , y = y ,color = Phase),size = 0.9 ) +
	scale_color_manual(values=colors)+
	labs(x="UMAP_1", y="UMAP_2") + 
	guides(color=guide_legend(override.aes = list(size=3))) +
	my_theme

ggsave(p,file=paste0(outpre,".pdf"),width = 8, height = 6)

