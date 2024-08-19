#-*-coding:utf-8-*-

def main(loom,outpre,color_map=None):
### 导入数据
	adata = scv.read(loom, cache=False) ## 读取loom文件

### 读绘图配色
	palette = None
	if color_map is not None:
		file = open(color_map, 'r')
		palette = file.read().split()
		file.close()

### 数据预处理
	scv.pp.filter_and_normalize(adata) ## 过滤和归一化
	scv.pp.moments(adata, n_pcs=30, n_neighbors=30) ## 速度估值的一阶和二阶距
### 速率分析
	scv.tl.velocity(adata) ## 速率分析
	scv.tl.velocity_graph(adata)  ##速度图数据
### 可视化
	scv.pl.velocity_embedding(adata)

### 绘图
	scv.pl.velocity_embedding_stream(adata,basis = 'umap',color = 'clusters',arrow_size=1.5,xlabel = "UMAP_1",ylabel = "UMAP_2", size = 35, alpha = 1, title = "",figsize = (12,10),legend_loc = 'right margin',save="%s.umap.velocity_stream.svg" % outpre, palette=palette)
	scv.pl.velocity_embedding_stream(adata,basis = 'umap',color = 'clusters',arrow_size=1.5,xlabel = "UMAP_1",ylabel = "UMAP_2", size = 35, alpha = 1, title = "",figsize = (12,10),legend_loc = 'on data',save="%s.umap.velocity_stream.on_data.svg" % outpre, palette=palette)


if __name__=="__main__":
	import sys,re,os

	if len(sys.argv) < 3:
		print("usage : python scVelo.py <loom> <outpre> (<color_map>)\n")
		sys.exit()

	import scvelo as scv
	import scanpy as sc
	import matplotlib
	import matplotlib.pyplot as plt
	plt.switch_backend('agg')
	matplotlib.use('Agg')

	main(*sys.argv[1:])
	print("done!")
