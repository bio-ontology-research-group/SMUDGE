
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation,Dropout
from sklearn.cross_validation import train_test_split
from sklearn.cross_validation import KFold
from keras.callbacks import ModelCheckpoint, EarlyStopping

import keras
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score

import pandas as pd
import itertools
import pdb
import json

import tensorflow as tf
import random


np.random.seed(33)
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
random.seed(12345)


data_dir = '../../Documents/smudge_data/'

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



gene_embeddings = {}
disease_embeddings = {}

data = pd.read_csv(data_dir+'P-Vec_embds_dfs_allclasses_32.txt', header = None, skiprows = 1, sep = ' ')
data = data.values
embds_dict = dict(zip(data[:,0],data[:,1:]))


for item in embds_dict:
	if item.isdigit(): 
	    gene_embeddings[item] = np.array(embds_dict[item], dtype = 'float32')
	elif item.startswith('OMIM:'):
		disease_embeddings[item] = np.array(embds_dict[item], dtype = 'float32')
	else:
	    continue

with open(data_dir+'disease_genes_human_hpo_sub.dict','r') as f:
	disease_genes = json.load(f)

disease_set = list(set(disease_genes.keys()).intersection(set(disease_embeddings.keys())))

x_train, x_test = train_test_split(disease_set, test_size=0.2, random_state=1)


dis_gene_tr = {}
dis_gene_ts = {}


for dis in x_train:
	if dis in disease_embeddings:
		dis_embds = disease_embeddings[dis]
		genes = disease_genes[dis]
		for gene in genes:
			if gene in gene_embeddings:
				gene_embds = gene_embeddings[gene]
				dis_gene_tr[(dis,gene)] = np.concatenate((dis_embds,gene_embds), axis=0)


for dis in x_test:
	if dis in disease_embeddings:
		dis_embds = disease_embeddings[dis]
		genes = disease_genes[dis]
		for gene in genes:
			if gene in gene_embeddings:
				gene_embds = gene_embeddings[gene]
				dis_gene_ts[(dis,gene)] = np.concatenate((dis_embds,gene_embds), axis=0)


diseases = disease_embeddings.keys()
genes = gene_embeddings.keys()
allpairs = list(itertools.product(set(diseases),set(genes)))
positive_pairs_tr = dis_gene_tr.keys()
positive_pairs_ts = dis_gene_ts.keys()


negative_pairs = set(allpairs) - set(positive_pairs_tr)
negative_pairs = set(negative_pairs) - set(positive_pairs_ts)

negative_pairs_tr = random.sample(negative_pairs, len(positive_pairs_tr))
negative_pairs = set(negative_pairs) - set(negative_pairs_tr)

negative_pairs_ts = random.sample(negative_pairs, len(positive_pairs_ts))


non_dis_gene_tr = {}
non_dis_gene_ts = {}

for (dis,gene) in negative_pairs_tr:
	dis_embds = disease_embeddings[dis]
	gene_embds = gene_embeddings[gene]
	non_dis_gene_tr[(dis,gene)] = np.concatenate((dis_embds,gene_embds), axis=0)


for (dis,gene) in negative_pairs_ts:
	dis_embds = disease_embeddings[dis]
	gene_embds = gene_embeddings[gene]
	non_dis_gene_ts[(dis,gene)] = np.concatenate((dis_embds,gene_embds), axis=0)


pos_tr = np.array(dis_gene_tr.values())
neg_tr = np.array(non_dis_gene_tr.values())

pos_ts = np.array(dis_gene_ts.values())
neg_ts = np.array(non_dis_gene_ts.values())

train_data = np.concatenate((pos_tr, neg_tr))
train_labels = np.concatenate((np.ones(len(dis_gene_tr), dtype = 'int32'),np.zeros(len(non_dis_gene_tr), dtype='int32')), axis=0)
tr_labels = keras.utils.to_categorical(train_labels, num_classes=None)

test_data = np.concatenate((pos_ts, neg_ts))
test_labels = np.concatenate((np.ones(len(dis_gene_ts), dtype = 'int32'),np.zeros(len(non_dis_gene_ts), dtype='int32')), axis=0)
ts_labels = keras.utils.to_categorical(test_labels, num_classes=None)

tf.set_random_seed(33)
model = Sequential()
model.add(Dense(2500, input_dim=64, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='sigmoid'))

model.compile(loss='binary_crossentropy',
          optimizer='rmsprop',
          metrics=['accuracy'])


earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
model.fit(train_data, tr_labels,validation_split = 0.1,epochs=100,callbacks=[earlystopper])

pred = model.predict_proba(test_data)[:,1]
test_labels = np.nonzero(ts_labels)[1]
print 'classifier prob thresholds: {}'.format(roc_auc_score(test_labels, pred))


dis_dict = {}
#apply the model for each disease
for dis in x_test:
	disemds = disease_embeddings[dis]
	if dis in disease_genes:
		gda = disease_genes[dis]
		gda = filter(None, gda)
		gda = list(set(gda))
		common  = set(genes).intersection(gda)

		test_emds = []
		for gene in genes:
			genemds = gene_embeddings[gene]
			pair_embds = np.concatenate((disemds, genemds), axis=0)
			test_emds.append(pair_embds)

		test_emds = np.array(test_emds, dtype='float32')
		y_pred = model.predict_proba(test_emds)[:,1]
		sorted_idx = np.argsort(y_pred)[::-1]
		sort_gene = [genes[arg] for arg in sorted_idx]

		label_vec = [0]*len(sort_gene)
		for gene in gda:
			if gene in sort_gene:
				label_vec[sort_gene.index(gene)] = 1

		dis_dict[dis] = label_vec
				
array_tp = np.zeros((len(dis_dict), len(genes)),dtype='float32')
array_fp = np.zeros((len(dis_dict), len(genes)), dtype = 'float32')

for i,row in enumerate(dis_dict.values()):
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
print('Number of Disease: {}'.format(len(dis_dict)))
print('Number of genes: {}'.format(len(genes)))
print('classifier rank thresholds: {}'.format(auc(fpr_r, tpr_r)))

pdb.set_trace()
