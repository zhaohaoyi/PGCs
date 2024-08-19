library(ggplot2)
library(reshape2)
library(scales)
library(dplyr)
library(RColorBrewer)
library(showtext)

library(ggcor)
# default theme settings
showtext_auto(enable = TRUE)

mytheme <- theme_cor() +
theme(
    text = element_text(family = "Arial"),

    axis.text  = element_text(color = "#000000", size = 12),

    legend.title = element_blank(),
    legend.text = element_text(hjust = 0.5, vjust = 0.5, margin = margin(2, 0, 0, 0)),

    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
)


readTable <- function (file, rowcol) {
    data <- read.table(file, head = T, sep = "\t", check.names = F, 
    row.names = 1, quote = "")
    if(rowcol == "row") {
        data <- t(data)
    }
    return(data)
}

file    <- "no_PGC_5D_05_DAZL.AllGene.avg_exp.xls"
outname <- "Cluster.cor.heatmap"

title  <- ""
xlab   <- ""
ylab   <- ""
islabel  <- as.logical("T")
islegend <- as.logical("T")

width  <- 7
height <- 7
colors <- c('#2F70AD','#FFFFFF','#BA2831')

# calculate correlation
data <- readTable(file, "col")
data.cor <- correlate(data, method = "pearson", cor.test = F)
data.cor.tb <- fortify_cor(data, method = "pearson", type = "full", cor.test = F)

write.table(data.frame("Sample" = rownames(data.cor$r), data.cor$r, check.names = F), 
    file = paste0(outdir, "/", outname, ".correlation.xls"), sep = "\t", quote = F, row.names = F)

width = 8
height = 8
axis.size = 10

# plot
p <- ggcor(data.cor.tb) +
     geom_square(aes(r0 = r, fill = r)) +
     scale_fill_gradientn(colors = colors) +
     coord_fixed(expand = FALSE) +
     geom_panel_grid(size = 1, color = "grey50") + 
     geom_panel_grid(size = 0.6, color = "#FFFFFF") + 
     labs(x = xlab, y = ylab, title = title) + 
     mytheme + 
     theme(axis.text = element_text(size = axis.size))

p <- p + geom_number(aes(num = r), size = text.size, color = "grey10", digits = 3)

ggsave(file = paste0(outname, ".pdf"), plot = p, width = width, height = height)

