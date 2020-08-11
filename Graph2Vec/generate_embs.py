import networkx as nx
import os
import karateclub as kc
import math
from numpy import asarray,savetxt

folder = '../../Networks/Processing/Networks'

files = os.listdir(folder)
files = list(filter(lambda a:a[:3]!='lfr',files))


dico_lab = {'fw':0, 'el':1,'st':2,'so':3}
lGraphs = []

sorted_name = ""
for file in files:
    print(file)
    G = nx.read_edgelist(folder+file,create_using=nx.DiGraph,nodetype = int)
    InDeg = G.in_degree()
    OutDeg = G.out_degree()
    lnodes = G.nodes()
    for node in lnodes:
        G.nodes[node]["feature"] = [InDeg[node], OutDeg[node]]
    G = nx.Graph(G)
    lGraphs.append(G)
    sorted_name += file[:-4]+'\n'
with open('sorted_name.txt','w') as f:
    f.write(sorted_name)


t_emb = [2**k for k in range(1,10)]
wls = range(1,5)
for wl in wls:
    for d_emb in t_emb:
        print(wl,'--',d_emb)
        model = kc.Graph2Vec(wl_iterations = wl, dimensions=d_emb, attributed=True,epochs = 100)
        model.fit(lGraphs)
        X = model.get_embedding()
        savetxt('Embeddings/wl_'+str(wl)+'_emb_'+str(d_emb)+'.csv',X,delimiter=",")
