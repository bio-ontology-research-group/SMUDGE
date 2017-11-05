import networkx as nx
import random
import pdb
import numpy as np
from gensim.models import Word2Vec



def read_into_dict(file_):
  some_dict = {}
  with open(file_) as f:
    for line in f:
      items = line.strip().split()
      item1 = items[0]
      item2 = items[1:]
      some_dict[item1] = item2

  return some_dict


def smudge_phenos(graph, num_walk, len_walk, genes_without_phenos, super_classes, genes_phenos, disease_phenos):
    
    all_walks = []
    nodes = graph.nodes()
    genes_with_phenos = genes_phenos.keys()
    genes_nodes = [node for node in nodes if node.isdigit()]
    disease_nodes = [node for node in nodes if node.startswith('OMIM:')]
    
    for dis in disease_nodes: #walk over disease nodes
      for walk in range(num_walk):
          path = []
          if dis in disease_phenos:
            phenos = disease_phenos[dis]

            pheno = random.choice(phenos)
            path2root = super_classes[pheno]
            path.append(dis)
            path.append(pheno)
            path.extend(path2root)
            all_walks.append(path)

    for node in genes_nodes: #walk over PPI nodes
        for walk in range(num_walk):
          path = []
          if node in genes_phenos:
              phenos = genes_phenos[node]
              pheno = random.choice(phenos)
              if pheno in super_classes:
                  path2root = super_classes[pheno]
                  path.append(node) #gene with pheno
                  path.append(pheno) #pheno
                  path.extend(path2root)
                  all_walks.append(path)
              
          elif node in genes_without_phenos: 
              neigbors = nx.single_source_shortest_path_length(graph ,source=node, cutoff=len_walk)
              genos_adj = list(set(neigbors).intersection(set(genes_with_phenos)))
              if not genos_adj:
                  continue   
              else:     
                  geno = random.choice(genos_adj)
                  phenos = genes_phenos[geno]
                  pheno = random.choice(phenos)
                  if pheno in super_classes:
                    path2root = super_classes[pheno]
                    path.append(node) #gene w/o pheno
                    path.append(geno) #gene with pheno
                    path.append(pheno) #pheno
                    path.extend(path2root)
                    all_walks.append(path)

    return all_walks



if __name__ == '__main__':

    data = '../../../Documents/smudge_data/'

    graph = nx.read_edgelist(data+'mouse_graph.txt', create_using=nx.DiGraph(), data=(('label', str),))
    genes_phenos = read_into_dict(data+'human_genes_mouse_phen.txt')
    disease_phenos = read_into_dict(data+'omim_hpos.txt')
    super_classes = read_into_dict(data+'phenomNet_super_classes.txt')
    genes_with_phenos = list(open(data+'genes_with_mouse_annotations.txt').readlines())
    genes_with_phenos = [item.strip() for item in genes_with_phenos]
    genes_without_phenos = list(open(data+'genes_without_mouse_phenos.txt').readlines())
    genes_without_phenos = [item.strip() for item in genes_without_phenos]

    print('The number of nodes in graph is: {}'.format(len(graph.nodes())))
    print('Walking PPIs and Phenotypes ...')

    walks = smudge_phenos(graph, 500, 20, genes_without_phenos, super_classes, genes_phenos, disease_phenos)
    print('Training the graph corpus...')
    model = Word2Vec(walks,size=512, window=10, min_count=1, sg =1, workers=24)

    model.save_word2vec_format(data+'mouse_embeddings_500_1.txt')

    pdb.set_trace()


