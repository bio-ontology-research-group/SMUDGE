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


with open(data_dir+'gwas_catalog_v1.0-associations_e91_r2018-01-16.tsv') as f:
	for line in f:
		if line.startswith('DATE'):
			continue

		items = line.strip().split('\t')
		disease = items[7]
		reported_genes = items[13].split(',')
		# efo = items[-2].split('/')[-1]
		genes_entr = []
		for item in reported_genes:
			if item in sym2id:
				genes_entr.append(sym2id[item])

		genes_entr = list(set(genes_entr))
		if disease in disease_genes:
			disease_genes[disease].extend(genes_entr)
		else:
			disease_genes[disease] = genes_entr

dis_gen_dict = {}
for item in disease_genes:
	genes = list(set(disease_genes[item]))
	dis_gen_dict[item] = genes

genes = gene_embeddings.keys()
genes_embeds = np.array(gene_embeddings.values())
gene_pairs = list(itertools.combinations(genes,2)) #all genes pairs

cosine_mat = cosine_similarity(genes_embeds)

#get upper triangle, exclude diagonal itself (self similarities)
cosine_pairs = cosine_mat[np.triu_indices(len(cosine_mat),1)]
cosine_idx = np.argsort(cosine_pairs)[::-1]
sort_pairs = [gene_pairs[arg] for arg in cosine_idx]
pairs_sorted_idx = dict(zip(sort_pairs,cosine_idx))

label_mat = dict()
for dis in dis_gen_dict:
	print dis
	dis_genes = dis_gen_dict[dis]
	print 'number of genes: {}'.format(len(dis_genes))
	spec_pairs = list(itertools.combinations(dis_genes,2)) #genes pairs of particular disease
	label_vec = [0]*len(sort_pairs)
	for pair in spec_pairs:
	 	if pair not in pairs_sorted_idx:
			pair1 = (pair[1],pair[0])
			if pair1 in pairs_sorted_idx:
				idx = pairs_sorted_idx[pair1]
				label_vec[idx] = 1
		else:
			idx = pairs_sorted_idx[pair]
			label_vec[idx] = 1

	label_mat[dis] = label_vec
	

array_tp = np.zeros((len(label_mat),len(gene_pairs)),dtype='float32')
array_fp = np.zeros((len(label_mat), len(gene_pairs)), dtype = 'float32')

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
print('Number of GWAS diseases: {} '.format(len(label_mat)))
print('Number of genes pairs ranked against: {}'.format(len(gene_pairs)))
print('auc {}'.format(auc(fpr_r, tpr_r)))

