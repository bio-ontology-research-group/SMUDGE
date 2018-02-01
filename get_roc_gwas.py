import pandas as pd
import pdb
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import auc
from numpy.linalg import norm
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json 



def negcum(rank_vec):
	rank_vec_cum = []
	prev = 0
	for x in rank_vec:
		if x == 0:
			x = x+1
			prev = prev + x
			rank_vec_cum.append(prev)
		else:
			rank_vec_cum.append(prev)
	rank_vec_cum = np.array(rank_vec_cum)
	return rank_vec_cum


data_dir = '../../Documents/smudge_data/'

data = pd.read_csv(data_dir+'human_embeddings_2000_hpo_512.txt', header = None, skiprows = 1, sep = ' ')
data = data.values
embds_dict = dict(zip(data[:,0],data[:,1:]))


gene_embeddings = {}
disease_embeddings = {}



for item in embds_dict:
	if item.isdigit(): 
		gene_embeddings[item] = np.array(embds_dict[item], dtype = 'float32')
	elif item.startswith('OMIM:'):
		disease_embeddings[item] = np.array(embds_dict[item], dtype = 'float32')
	else:
		continue

sym2id = {}
disease_genes = {}
with open(data_dir+'geneID_symbol-mod.txt') as f:
	for line in f:
		items = line.split()
		geneID = items[0]
		sym = items[1]
		sym2id[sym] = geneID


with open(data_dir+'gwas_catalog_v1.0.1-associations_e91_r2018-01-23.tsv') as f:
	for line in f:
		if line.startswith('DATE'):
			continue

		items = line.strip().split('\t')
		# disease_name = items[7]
		reported_genes = items[13].split(',')
		efo = items[-2].split('/')[-1]
		genes_entr = []
		for item in reported_genes:
			if item in sym2id:
				genes_entr.append(sym2id[item])

		genes_entr = list(set(genes_entr))
		if efo in disease_genes:
			disease_genes[efo].extend(genes_entr)
		else:
			disease_genes[efo] = genes_entr


genes = gene_embeddings.keys()
genes_embeds = np.array(gene_embeddings.values())


label_mat = dict()
for dis in disease_genes:
	dis_genes = disease_genes[dis]
	common_genes = list(set(genes).intersection(set(dis_genes)))

	if len(common_genes) <= 1:
		continue
	
	print dis
	print 'number of genes: {}'.format(len(common_genes))
	gen_dict = {}
	for gene1 in common_genes:
		gene_embds = gene_embeddings[gene1]
		gen_embds = list()
		gen_embds.append(gene_embds)
		cos_sim = cosine_similarity(gen_embds, genes_embeds)
		cos_sim = cos_sim.flatten()
		sort_sim = np.argsort(cos_sim)[::-1]
		sort_genes = [genes[arg] for arg in sort_sim]
		label_vec = [0]*len(genes)
		for gene2 in common_genes:
			if gene1 != gene2:
				if gene2 in sort_genes:
						label_vec[sort_genes.index(gene2)] = 1

		gen_dict[gene1] = label_vec
		
	labels = np.array(gen_dict.values())
	labels = np.max(labels, axis=0) #get all positives pairs labels for specific dis, all others negative
	label_mat[dis] = labels
	

array_tp = np.zeros((len(label_mat), len(genes)),dtype='float32')
array_fp = np.zeros((len(label_mat), len(genes)), dtype = 'float32')

for i,row in enumerate(label_mat.values()):
	elem = np.asarray(row, dtype='float32')
	tpcum = np.cumsum(elem)
	fpcum = negcum(elem)  	
	array_tp[i] = tpcum
	array_fp[i] = fpcum	

#compute fpr and tpr Rob's way 
tpsum = np.sum(array_tp, axis = 0)
fpsum = np.sum(array_fp, axis = 0)
tpr_r = tpsum/max(tpsum)
fpr_r = fpsum/max(fpsum)
auc_data2 = np.c_[fpr_r, tpr_r]
file1 = open(data_dir+'gwas_results.txt','w')

file1.write('Number of GWAS diseases: {} \n'.format(len(label_mat)))
file1.write('Number of genes ranked against: {}\n'.format(len(genes)))
file1.write('auc {}\n'.format(auc(fpr_r, tpr_r)))
file1.close()
pdb.set_trace()

