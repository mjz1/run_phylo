import json
import pandas as pd
import os 
import argparse
import subprocess
import glob
import re
import gzip
import zipfile
from functools import reduce
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri

def convert_rda(outdir, paths):
#	pandas2ri.activate()
	rdas = []
	for path in paths: 
		#### TEMP ####
		temp_path = os.path.join(outdir, os.path.splitext(os.path.basename(path))[0]+".tmp") 
		command = "Rscript /hpf/largeprojects/adam/projects/lfs/run_phylowgs/rda_to_tab.R --mutect " + path + " --outfile " + temp_path
		subprocess.call(command, shell=True)
		df = pd.read_csv(temp_path, header=0, delimiter='\t')
		##############
	#	rda = robjects.r.load(path)[0]
	#	df = robjects.r[rda]
		rdas.append(df)
	return(rdas)


## use summ to calculate best tree 
def get_best_tree(outdir, summ_path):
	with open(summ_path) as f:
		summ_dict = json.load(f)
	densities = pd.DataFrame(summ_dict['tree_densities'], index=[0]).transpose()
	max_density = densities.max()[0]
	best_tidx = densities.idxmax()[0]
	print("The best tree is tree {} with a density of {}".format(best_tidx, max_density))
	best_tree = glob.glob(os.path.join(outdir, "results","mutass", best_tidx+".json"))[0]
	return(best_tree)

## parse input
def parse_input(outdir):
	if not os.path.isdir(os.path.join(outdir ,"results", "mutass")):
		zip_ref = zipfile.ZipFile(glob.glob(os.path.join(outdir ,"results/*.mutass.zip"))[0], 'r')
		zip_ref.extractall(os.path.join(outdir ,"results","mutass"))
		zip_ref.close()
	if glob.glob(os.path.join(outdir ,"results/*.muts.json.gz")):
		gunzip_muts = "gunzip " + glob.glob(os.path.join(outdir ,"results/*.muts.json.gz"))[0]
		subprocess.call(gunzip_muts, shell=True)
	if glob.glob(os.path.join(outdir ,"results/*.summ.json.gz")):
		gunzip_summs = "gunzip " + glob.glob(os.path.join(outdir, "results/*.summ.json.gz"))[0]
		subprocess.call(gunzip_summs, shell=True)
	path_names = ['cnv_physical_path', 'ssm_physical_path', 'annotated_ssm_path', 'logical_ids_path', 'summ_path']
	paths = [glob.glob(os.path.join(outdir,'*_cnv_data.txt'))[0], glob.glob(os.path.join(outdir,'*_ssm_data.txt'))[0], os.path.join(outdir , 'tmp/'), glob.glob(os.path.join(outdir ,"results/*.muts.json"))[0], glob.glob(os.path.join(outdir, "results/*.summ.json"))[0]]	
	paths_dict = dict(zip(path_names, paths))
	return(paths_dict)

## create table for best tree
def create_table(outdir, mutect2):

	paths = parse_input(outdir)
	best_tree_path = get_best_tree(outdir, paths['summ_path'])

	with open(best_tree_path) as f:
	    phylo_json_dict = json.load(f)
	samples = phylo_json_dict['dataset_name'].split(',')
	phylo_muts = pd.DataFrame(phylo_json_dict['mut_assignments']).transpose()
	muts = pd.concat([pd.DataFrame(phylo_muts.cnvs.values.tolist(),index= phylo_muts.index), pd.DataFrame(phylo_muts.ssms.values.tolist(), index= phylo_muts.index)], axis=1)
	muts = pd.melt(muts.transpose()).dropna()
	phylo_dict = dict(zip(muts["value"], muts["variable"]))


	cnv_physical = pd.read_csv(paths['cnv_physical_path'], header=0, delimiter='\t')
	cnv_physical[['chrom','start','end', 'major_cn','minor_cn','cell_prev']] = cnv_physical['physical_cnvs'].str.split(',',expand=True)
	cnv_physical['gene'] = cnv_physical['physical_cnvs'].str.replace('chrom=|start=|end=|major_cn=|minor_cn=|cell_prev=','').str.split(',', expand=True).loc[: ,0:2].apply(lambda x:'%s_%s_%s' % (x[0],x[1],x[2]),axis=1)
	cnv_physical[[s + "_a" for s in samples]] = cnv_physical['a'].str.split(',', expand=True)
	cnv_physical[[s + "_d" for s in samples]] = cnv_physical['d'].str.split(',', expand=True)
	ssm_physical = pd.read_csv(paths['ssm_physical_path'], header=0, delimiter='\t')
	ssm_physical[[s + "_a" for s in samples]] = ssm_physical['a'].str.split(',', expand=True)
	ssm_physical[[s + "_d" for s in samples]] = ssm_physical['d'].str.split(',', expand=True)	
	#ssm_dict = ssm_physical[['id'] + [s + "_a" for s in samples]].groupby('id').apply(lambda dfg: dfg.to_dict(orient='list')).to_dict()
	variants = pd.concat([ssm_physical[['id', 'gene']], cnv_physical[['cnv', 'gene']].rename(columns={'cnv':'id'})], ignore_index=True)
	variants['clone'] = variants['id'].map(phylo_dict)

	annotated_ssms = convert_rda(outdir, mutect2)

	for annotated_ssm in annotated_ssms:
		annotated_ssm['id'] = annotated_ssm['snvid'].str.split('_', expand=True).loc[: ,0:1].apply(lambda x:'%s_%s' % (x[0],x[1]),axis=1)
		annotated_ssm_dict = dict(zip(annotated_ssm["id"],annotated_ssm["snvid"]))
		variants['gene'] = variants['gene'].replace(annotated_ssm_dict)

	## a = ref_reads ; d = total_reads
	variants_a = pd.concat([ssm_physical[['id'] + [s + "_a" for s in samples]], cnv_physical[['cnv'] + [s + "_a" for s in samples]].rename(columns={'cnv':'id'})], ignore_index=True)
	variants_d = pd.concat([ssm_physical[['id'] + [s + "_d" for s in samples]], cnv_physical[['cnv'] + [s + "_d" for s in samples]].rename(columns={'cnv':'id'})], ignore_index=True)
	dfs = [variants, variants_a, variants_d]
	variants = reduce(lambda left,right: pd.merge(left,right,on='id'), dfs)

	variants.rename(columns={'gene':'snvid'},inplace=True)
	
	for path, annotated_ssm in zip(mutect2, annotated_ssms):
		annotated_sample_path = os.path.join(outdir, "results", os.path.basename(path).split('_annotated_filtered_clipped.rda')[0]+"_phylowgs.txt")
		sample_variants = pd.merge(annotated_ssm, variants, how='left', on='snvid')
		sample_variants.to_csv(annotated_sample_path, index=False, sep='\t')

	annovar = list(filter(re.compile("annovar.*").match, annotated_ssms[0].columns))
	cols = annovar + ["ensembl_gene","hgnc_gene","cosmic_census","aa", "snvid"] 
	ssms_annotated = pd.concat(annotated_ssms)[cols].drop_duplicates(subset=['snvid'], keep='first')

	variants_multisample = pd.merge(variants, ssms_annotated, how='left', on='snvid')
	outfile = os.path.join(outdir, "results", samples[0].split('_')[0] + "_phylowgs.txt")
	variants_multisample.to_csv(outfile, index=False, sep='\t')
#	print(variants.groupby(['clone']).agg({'position':['count',lambda x: ','.join(x)],'id':[lambda y: y.str.contains('c').count(), lambda x: x.str.contains('s').count()]}))

def main():

	parser = argparse.ArgumentParser(
		description='Create formatted phylowgs output.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument('-o','--outdir', dest='outdir', default=os.getcwd(),
		help='Output destination for variants')
	parser.add_argument('-m','--mutect', nargs='+', dest='mutect2', help='Annotated rda for EACH sample in multisample', required=True)

	args = parser.parse_args()

	create_table(args.outdir, args.mutect2)
	

main()
