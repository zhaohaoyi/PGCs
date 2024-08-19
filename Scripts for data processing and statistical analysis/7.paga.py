if __name__ == "__main__":
	import sys
	import os
	import argparse
	parser = argparse.ArgumentParser(description='PAGA pipeline')
	parser.add_argument('--config', dest='config', required=True, help='parameters file in yaml format')
	parser.add_argument('--loom', dest='loom', required=True, help='seurat.loom file')
	parser.add_argument('--outdir', dest='outdir', required=True, help='output directory')
	args = parser.parse_args()

	libpath = os.path.dirname(__file__) + '/libs/python'
	sys.path.append(libpath)
	from scanpy_library import *
	sc.settings.verbosity = 3
	conf_file = args.config
	loom_file = args.loom
	outdir = args.outdir
	results_file = 'paga.h5ad'
	para = read_conf(conf_file)
	setwd(outdir)

	adata = seurat2adata(loom_file)
	
	adata = RunPAGA(adata, root_cluster=para['root_cluster'], root_cell = para['root_cell'])
	write_pd_table(
			adata.obs,
			output="Trajectory.Data.xls",
			first_colname="Cells",
			columns={'cluster':'Cluster', 'sample':'Sample','dpt_pseudotime':'Pseudotime'}
			)
	adata.write(results_file, compression='gzip')


