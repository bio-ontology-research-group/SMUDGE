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


def smudge_phenos(graph, num_walk, genes_without_phenos, super_classes, genes_phenos, disease_phenos):
    
    all_walks = []
    nodes = graph.nodes()
    genes_with_phenos = genes_phenos.keys()
    genes_nodes = [node for node in nodes if node.isdigit()]
    disease_nodes = [node for node in nodes if node.startswith('OMIM:')]
    print('Number of genes nodes: {}'.format(len(genes_nodes)))
    print('Number of disease nodes: {}'.format(len(disease_nodes)))


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
          neigbors = graph.neighbors(node)
          ppi_adj = list(set(neigbors).intersection(genes_nodes))
          if ppi_adj:
              rand_node = random.choice(ppi_adj)
              stumbles = 0

              while(rand_node not in genes_with_phenos): #keep walking until find gene/w pheno
                  if (stumbles > 5):
                      break

                  neigbors = graph.neighbors(rand_node) 
                  if neigbors:
                      rand_node = random.choice(neigbors)
                      stumbles = stumbles + 1

              if stumbles > 5:
                  continue  

              phenos = genes_phenos[rand_node]
              pheno = random.choice(phenos)
              if pheno in super_classes:
                path2root = super_classes[pheno]
                path.append(node)
                path.append(pheno) #pheno
                path.extend(path2root)
                all_walks.append(path)


    return all_walks



if __name__ == '__main__':

    data = '../../../Documents/smudge_data/'

    graph = nx.read_edgelist(data+'human_graph.txt', create_using=nx.DiGraph(), data=(('label', str),))
    genes_phenos = read_into_dict(data+'human_genes_hpos.txt')
    disease_phenos = read_into_dict(data+'omim_hpos.txt')
    super_classes = read_into_dict(data+'hp_super_classes.txt')
    genes_without_phenos = list(open(data+'genes_without_human_phenos.txt').readlines())
    genes_without_phenos = [item.strip() for item in genes_without_phenos]

    print('The number of nodes in graph is: {}'.format(len(graph.nodes())))
    print('Walking PPIs and Phenotypes ...')

    walks = smudge_phenos(graph, 500, genes_without_phenos, super_classes, genes_phenos, disease_phenos)
    print('Training the graph corpus...')
    model = Word2Vec(walks,size=512, window=10, min_count=1, sg =1, workers=24)

    model.save_word2vec_format(data+'human_embeddings_500.txt')

    pdb.set_trace()


