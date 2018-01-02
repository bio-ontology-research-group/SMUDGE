
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation,Dropout
from sklearn.cross_validation import train_test_split
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
# one thread to ensure reproducibility 
session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
random.seed(12345)


data_dir = '../../Documents/smudge_data/'
data = pd.read_csv(data_dir+'human_embeddings_500_hpo.txt', header = None, skiprows = 1, sep = ' ')
data = data.values


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

for item in data:
	it = item[0].split('/')[-1]
   	if it.isdigit(): 
   	    gene_embeddings[it] = np.array(item[1:], dtype = 'float32')
   	elif it.startswith('OMIM:'):
   	    disease_embeddings[it] = np.array(item[1:], dtype = 'float32')
   	else:
   	    continue


with open(data_dir+'disease_genes_human_hpo.dict','r') as f:
	disease_genes = json.load(f)

dis_gene_embds = {}
dis_gene = {}

for dis in disease_genes:
	if dis in disease_embeddings:
		dis_embds = disease_embeddings[dis]
		genes = disease_genes[dis]
		for gene in genes:
			if gene in gene_embeddings:
				gene_embds = gene_embeddings[gene]
				dis_gene_embds[(dis,gene)] = np.concatenate((dis_embds,gene_embds), axis=0)


diseases = disease_embeddings.keys()
genes = gene_embeddings.keys()
allpairs = list(itertools.product(set(diseases),set(genes)))
positive_pairs = dis_gene_embds.keys()

negative_pairs = set(allpairs) - set(positive_pairs)

negative_pairs_ = random.sample(negative_pairs, len(positive_pairs))
non_assoc_embds = {}

for item in negative_pairs_:
	dis_embds = disease_embeddings[item[0]]
	gene_embds = gene_embeddings[item[1]]
	non_assoc_embds[(item[0],item[1])] = np.concatenate((dis_embds,gene_embds), axis=0)


positive_pairs = np.array(dis_gene_embds.keys())
negative_pairs = np.array(non_assoc_embds.keys())
data_pairs = np.concatenate((positive_pairs, negative_pairs))
labels = np.concatenate((np.ones(len(positive_pairs), dtype = 'int32'),np.zeros(len(negative_pairs), dtype='int32')), axis=0)

one_hot_labels = keras.utils.to_categorical(labels, num_classes=None)
x_train, x_test, y_train, y_test = train_test_split(data_pairs, one_hot_labels, test_size=0.2, random_state=1)

train_data = []
test_data = []
for item in x_train:
	disemd = disease_embeddings[item[0]]
	genembd = gene_embeddings[item[1]]
	pair_embds = np.concatenate((disemd, genembd), axis=0)
	train_data.append(pair_embds)

pos_test = []
for item in x_test:
	disemd = disease_embeddings[item[0]]
	genembd = gene_embeddings[item[1]]
	pair_embds = np.concatenate((disemd, genembd), axis=0)
	test_data.append(pair_embds)	
	if item in positive_pairs:
		pos_test.append(item)

train_data = np.array(train_data, dtype = 'float32')
test_data = np.array(test_data, dtype = 'float32')

# for optimizing ann hidden units,
#file1 = open('hidden_units.txt','w')
# for ii in range(100,3100,100):

tf.set_random_seed(33)
model = Sequential()
model.add(Dense(2000, input_dim=1024, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(2, activation='sigmoid'))

model.compile(loss='binary_crossentropy',
          optimizer='rmsprop',
          metrics=['accuracy'])


earlystopper = EarlyStopping(monitor='val_loss', patience=5, verbose=1)
model.fit(train_data, y_train,validation_split = 0.2,epochs=100,callbacks=[earlystopper])


label_mat = {}
test_pair = []
#apply the model for each disease, exclude genes in train
for dis in pos_test:
	disemds = disease_embeddings[dis[0]]
	if dis[0] in disease_genes:
		gda = disease_genes[dis[0]]
		gda = filter(None, gda)
		gda = list(set(gda))
		common  = set(genes).intersection(gda)

		if len(common) > 1: 	
			train_genes = [item[1] for item in x_train if item[0] == dis[0]]				
			to_retrieve = set(common) - set(train_genes)
			to_test = list(set(genes) - set(train_genes))
			test_emds = []
			if len(to_retrieve) >= 1:
				for gene in to_test:
					genemds = gene_embeddings[gene]
					pair_embds = np.concatenate((disemds, genemds), axis=0)
					test_emds.append(pair_embds)

				for gen in to_retrieve:
					test_pair.append([dis[0],gen])

				test_emds = np.array(test_emds, dtype='float32')
				y_pred = model.predict_proba(test_emds)[:,1]
				sorted_idx = np.argsort(y_pred)[::-1]
				sort_gene = [to_test[arg] for arg in sorted_idx]

				label_vec = [0]*len(sort_gene)
				for gene in to_retrieve:
					if gene in sort_gene:
						label_vec[sort_gene.index(gene)] = 1
				label_mat[dis[0]] = label_vec
				

#get max label vec dimension to make all disease ranks dimesnion is same
#every disease will have differnet numbers of genes to test with after removing
# genes in training set

# for all diseases
len_vec = []
for item in label_mat:
	len_vec.append(len(label_mat[item]))

col_num = max(len_vec)
array_tp = np.zeros((len(label_mat), col_num),dtype='float32')
array_fp = np.zeros((len(label_mat), col_num), dtype = 'float32')

for i,row in enumerate(label_mat.values()):
        elem = np.asarray(row, dtype='float32')
        tofill = col_num - len(row)
        tpcum = np.cumsum(elem)
        tpcum = np.append(tpcum,np.ones(tofill)*tpcum[-1])
        fpcum = negcum(elem)
    	fpcum = np.append(fpcum,np.ones(tofill)*fpcum[-1])
        array_tp[i] = tpcum
        array_fp[i] = fpcum


#compute fpr and tpr Rob's way 
tpsum = np.sum(array_tp, axis = 0)
fpsum = np.sum(array_fp, axis = 0)
tpr_r = tpsum/max(tpsum)
fpr_r = fpsum/max(fpsum)
print('Number of Disease: {}'.format(len(label_mat)))
print('Number of genes: {}'.format(len(genes)))
print('classifier rank thresholds: {}'.format(auc(fpr_r, tpr_r)))
# file1.write('{} units: {}\n'.format(ii,auc(fpr_r, tpr_r)))
# file1.close()

auc_data = np.c_[fpr_r, tpr_r]
np.savetxt('smudge_results/classifer_rank_ann_smudge'+'.txt', auc_data, fmt = '%s')


pdb.set_trace()
