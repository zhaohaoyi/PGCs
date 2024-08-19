library(pheatmap)

infile = "target.gene.heatmap.xls"
outpfx = "target.gene.heatmap"

data = read.table(infile,header=T,row.names=1,sep="\t",check.names = F,quote = "",comment.char = "")

## deal color
incolors = "#2F70AD,#FFFFFF,#BA2831"
color=(unlist(strsplit(incolors,',')))
mycolors=colorRampPalette(color)(100)

# plot
outfile=paste(outpfx,".pdf",sep="")
pheatmap(data,
	scale="row",
	cluster_rows = FALSE,
	cluster_cols = FALSE,
	color = mycolors,
	cellwidth = 10,
	fontsize = 13,
	show_rownames = TRUE,
	show_colnames = TRUE,
	filename = outfile
)

