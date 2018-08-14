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
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

def parse_phylotab(tab_path, parentdir):
	tab = pd.read_csv(tab_path, header=0, delimiter='\t')
	grouped = tab[['SNV_File','Multi_sample']].groupby('Multi_sample', as_index=False).agg({'SNV_File':[lambda x: ','.join(str(y) for y in x)]})
	for sample in grouped['Multi_sample']:
		sampledir = glob.glob(os.path.join(parentdir, sample+"_phylowgs.multi"))
		if sampledir:
			print("Creating phylowgs output for {}".format(sample))
			mutect = grouped.loc[grouped.Multi_sample == sample, 'SNV_File'].values[0][0].split(',')
			try:
				create_table(sampledir[0],mutect)
			except:
				continue
		else:
			print("ERROR: {} could not be found, skipping {}".format(sampledir, sample))
	
def convert_rda(outdir, paths):
	pandas2ri.activate()
	rdas = []
	for path in paths: 
		rda = robjects.r.load(path)[0]
		df = robjects.r[rda]
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
	if not os.path.isdir(os.path.join(outdir ,"results")):
		print("ERROR: {} not found, skipping sample...".format(os.path.join(outdir ,"results")))

	if glob.glob(os.path.join(outdir ,"results/*.muts.json.gz")):
		gunzip_muts = "gunzip " + glob.glob(os.path.join(outdir ,"results/*.muts.json.gz"))[0]
		subprocess.call(gunzip_muts, shell=True)
	elif glob.glob(os.path.join(outdir ,"results/*.muts.json")):
		pass
	else:
		print("ERROR: {} not found, skipping sample...".format(os.path.join(outdir ,"results/*.muts.json.gz")))

	if glob.glob(os.path.join(outdir ,"results/*.summ.json.gz")):
		gunzip_summs = "gunzip " + glob.glob(os.path.join(outdir, "results/*.summ.json.gz"))[0]
		subprocess.call(gunzip_summs, shell=True)
	elif glob.glob(os.path.join(outdir ,"results/*.summ.json")):
		pass
	else:
		print("ERROR: {} not found, skipping sample...".format(os.path.join(outdir ,"results/*.summ.json.gz")))

	if not os.path.isdir(os.path.join(outdir ,"results", "mutass")):
		zip_ref = zipfile.ZipFile(glob.glob(os.path.join(outdir ,"results/*.mutass.zip"))[0], 'r')
		zip_ref.extractall(os.path.join(outdir ,"results","mutass"))
		zip_ref.close()

	path_names = ['cnv_physical_path', 'ssm_physical_path', 'annotated_ssm_path', 'logical_ids_path', 'summ_path']
	paths = [glob.glob(os.path.join(outdir,'*_cnv_data.txt'))[0], glob.glob(os.path.join(outdir,'*_ssm_data.txt'))[0], os.path.join(outdir , 'tmp/'), glob.glob(os.path.join(outdir ,"results/*.muts.json"))[0], glob.glob(os.path.join(outdir, "results/*.summ.json"))[0]]	
	paths_dict = dict(zip(path_names, paths))
	return(paths_dict)

## create table for best tree
def create_table(outdir, mutect2):

	paths = parse_input(outdir)
	best_tree_path = get_best_tree(outdir, paths['summ_path'])

	best_tree(outdir)

	with open(best_tree_path) as f:
	    phylo_json_dict = json.load(f)
	
	with open(paths['summ_path']) as f:
            summ_json_dict = json.load(f)
	
	samples = summ_json_dict['params']['samples']
	phylo_muts = pd.DataFrame(phylo_json_dict['mut_assignments']).transpose()
	muts = pd.concat([pd.DataFrame(phylo_muts.cnvs.values.tolist(),index= phylo_muts.index), pd.DataFrame(phylo_muts.ssms.values.tolist(), index= phylo_muts.index)], axis=1)
	muts = pd.melt(muts.transpose()).dropna()
	phylo_dict = dict(zip(muts["value"], muts["variable"]))

	cnv_physical = pd.read_csv(paths['cnv_physical_path'], header=0, delimiter='\t', dtype={'a':str, 'd':str})
	physical_id = pd.concat([cnv_physical['physical_cnvs'].str.split(';',expand=True),cnv_physical[['cnv','a','d']]],axis=1); del(cnv_physical)
	physical_id = pd.melt(physical_id, id_vars=['cnv','a', 'd']).dropna()
	physical_id[['chrom','start','end', 'major_cn','minor_cn','cell_prev']] = physical_id['value'].str.split(',',expand=True)
	physical_id['value'] = physical_id['value'].str.replace('chrom=|start=|end=|major_cn=|minor_cn=|cell_prev=','').str.split(',', expand=True).loc[: ,0:2].apply(lambda x:'%s_%s_%s' % (x[0],x[1],x[2]),axis=1)
	physical_id[[s + "_a" for s in samples]] = physical_id['a'].str.split(',', expand=True)
	physical_id[[s + "_d" for s in samples]] = physical_id['d'].str.split(',', expand=True)
	ssm_physical = pd.read_csv(paths['ssm_physical_path'], header=0, delimiter='\t', dtype={'a':str, 'd':str})
	ssm_physical[[s + "_a" for s in samples]] = ssm_physical['a'].str.split(',', expand=True)
	ssm_physical[[s + "_d" for s in samples]] = ssm_physical['d'].str.split(',', expand=True)	
	
	#ssm_dict = ssm_physical[['id'] + [s + "_a" for s in samples]].groupby('id').apply(lambda dfg: dfg.to_dict(orient='list')).to_dict()
	variants = pd.concat([ssm_physical[['id', 'gene']], physical_id[['cnv', 'value']].rename(columns={'cnv':'id', 'value':'gene'})], ignore_index=True)
	variants['clone'] = variants['id'].map(phylo_dict)
	annotated_ssms = convert_rda(outdir, mutect2)

	for annotated_ssm in annotated_ssms:
		annotated_ssm['id'] = annotated_ssm['snvid'].str.split('_', expand=True).loc[: ,0:1].apply(lambda x:'%s_%s' % (x[0],x[1]),axis=1)
		annotated_ssm_dict = dict(zip(annotated_ssm["id"],annotated_ssm["snvid"]))
		variants['gene'] = variants['gene'].replace(annotated_ssm_dict)

	## a = ref_reads ; d = total_reads
	variants_a = pd.concat([ssm_physical[['id'] + [s + "_a" for s in samples]], physical_id[['cnv'] + [s + "_a" for s in samples]].rename(columns={'cnv':'id'})], ignore_index=True)
	variants_d = pd.concat([ssm_physical[['id'] + [s + "_d" for s in samples]], physical_id[['cnv'] + [s + "_d" for s in samples]].rename(columns={'cnv':'id'})], ignore_index=True)
	dfs = [variants, variants_a, variants_d]
	variants_final = reduce(lambda left,right: pd.merge(left,right,on='id'), dfs)
	variants_final.rename(columns={'gene':'snvid'},inplace=True)
	variants_final = variants_final.drop_duplicates()
	for path, annotated_ssm in zip(mutect2, annotated_ssms):
		annotated_sample_path = os.path.join(outdir, "results", os.path.basename(path).split('_annotated_filtered_clipped.rda')[0]+"_phylowgs.txt")
		sample_variants = pd.merge(annotated_ssm, variants_final, how='left', on='snvid')
		sample_variants.to_csv(annotated_sample_path, index=False, sep='\t')

	annovar = list(filter(re.compile("annovar.*").match, annotated_ssms[0].columns))
	cols = annovar + ["ensembl_gene","hgnc_gene","cosmic_census","aa", "snvid"] 
	ssms_annotated = pd.concat(annotated_ssms)[cols].drop_duplicates(subset=['snvid'], keep='first')
	variants_multisample = pd.merge(variants_final, ssms_annotated, how='left', on='snvid').drop_duplicates()
	outfile = os.path.join(outdir, "results", samples[0].split('_')[0] + "_phylowgs.txt")
	variants_multisample.to_csv(outfile, index=False, sep='\t')
#	print(variants.groupby(['clone']).agg({'position':['count',lambda x: ','.join(x)],'id':[lambda y: y.str.contains('c').count(), lambda x: x.str.contains('s').count()]}))

def best_tree(outdir):
	paths = parse_input(outdir)
	with open(paths['summ_path']) as f:
		summ_dict = json.load(f)
	densities = pd.DataFrame(summ_dict['tree_densities'], index=[0]).transpose()
	best_tidx = densities.idxmax()[0]
#	best_tidx="10"
	samples = summ_dict['params']['samples']
	best_tree = pd.DataFrame(summ_dict['trees'])[best_tidx]
	best_tree.loc[~best_tree.index.isin(['populations'])].to_csv(os.path.join(outdir, "best_tree.txt"), sep='\t')
	structure = best_tree['structure']
	print_best_tree(outdir, structure)
	populations = pd.DataFrame(best_tree['populations']).transpose()
	population = pd.concat({'Cellular Prevalence': populations.cellular_prevalence.apply(pd.Series), 'Variants':populations.iloc[: ,1:3]}, axis=1)
	population.columns.set_levels(samples +['num_cnvs','num_ssms'], level=1,inplace=True)
	population.insert(0, 'Population', population.index)
	population.to_csv(os.path.join(outdir, 'best_tree_populations.txt'),sep='\t', index=False)


def print_best_tree(outdir, structure):
	if not os.path.exists(os.path.join(outdir, "best_tree.txt")):
		fw = open(os.path.join(outdir, "best_tree.txt"), "w")
		levels = [value for values in structure.keys() for value in values]
		for level, nodes in structure.items():
			fw.write(' ----- ')
			for node in nodes:
				fw.write('%s' % node) 

def main():

	parser = argparse.ArgumentParser(
		description='Create formatted phylowgs output.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	
	parser.add_argument('-s','--sampledir', required=False, dest='sampledir', help='Directory containing sample phylowgs calls')
	parser.add_argument('-m','--mutect', required=False, nargs='+', dest='mutect2', help='Annotated rda for EACH sample in multisample')
	parser.add_argument('-t','--tab', required=False, dest='tab', help='Output destination for variants')
	parser.add_argument('-p','--parentdir', required=False, dest='parentdir', help='Parent directory containing ALL phylowgs calls in tab file')
	args = parser.parse_args()

	if args.tab and args.parentdir:
		parse_phylotab(args.tab, args.parentdir)
	elif args.sampledir and args.mutect2:
		create_table(args.sampledir, args.mutect2)
	else:
		parser.error('--outdir and --mutect must be given together OR --tab must be given')

main()

