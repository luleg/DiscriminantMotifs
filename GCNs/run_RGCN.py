import argparse
import time
import numpy as np
import networkx as nx
import torch
import torch.nn as nn
import torch.nn.functional as F
import dgl
from dgl.data import register_data_args
from dgl.nn.pytorch import GraphConv, HeteroGraphConv
from dgl.nn.pytorch.glob import SumPooling, AvgPooling, MaxPooling
import pickle as pkl
import argparse as arg

nbItes = 5
earlyStop=20
nbEpochs = 100
nbFeats = 2
nbClasses = 4
drop = 0.5
lr = 1e-2
weight_decay = 5e-4



parser = arg.ArgumentParser(usage="Run RGCN on data provided using create_graphtest.py. Write scores in a dict and save it as a pickle file. \n-pooling: Type of pooling layer to use after the GCLs (default:mean).\n-hidd :  Size of the embeddings between GCLs (default:12).\n-nbLay : Number of GCLs (default:2).",add_help=True)
parser.add_argument('-pooling',choices=['sum','mean','max'],default='mean')
parser.add_argument('-hidd',default='8')
parser.add_argument('-nbLay',default='2')
args = parser.parse_args()



nbLayers = int(args.nbLay)
nbHidden = int(args.hidd)
mpool=args.pooling

file2Save = 'RGCN_pool'+mpool+'_nbLay'+str(nbLayers)+'_hidd'+str(nbHidden)

dicoScores = {'tp':[],'fp':[],'fn':[],'stop':[],'time':[]}


seed = 123
np.random.seed(seed)
torch.manual_seed(seed)
##################################################################################

class Classifier(nn.Module):
    def __init__(self,
                 g,
                 rel_names,
                 in_feats,
                 n_hidden,
                 n_classes,
                 n_layers,
                 type_aggreg,
                 activation,
                 pooling,
                 dropout):
        super(Classifier, self).__init__()
        self.g = g
        self.activation = activation
        self.layers = nn.ModuleList()
        # input layer
        self.layers.append(HeteroGraphConv({
            rel:GraphConv(in_feats, n_hidden,allow_zero_in_degree=True,norm='both')
            for rel in rel_names
            },
            aggregate=type_aggreg
            ))
        # hidden layers
        for i in range(n_layers - 1):
            self.layers.append(HeteroGraphConv({
                rel:GraphConv(n_hidden, n_hidden, allow_zero_in_degree=True,norm='both')
                for rel in rel_names
                },
                aggregate=type_aggreg
                ))
        # output layer
        self.dropout = nn.Dropout(p=dropout)


        if pooling == 'sum':
            self.pool = SumPooling()
        elif pooling == 'mean':
            self.pool = AvgPooling()
        elif pooling == 'max':
            self.pool = MaxPooling()
        else:
            raise NotImplementedError

        self.classify = nn.Linear(n_hidden,n_classes)



    def forward(self, features):
        nbLay = len(self.layers)
        h = features
        for i, layer in enumerate(self.layers):
            if i != 0:
                h = {k:self.dropout(v) for k,v in h.items()}
            h = layer(self.g, h)
            h = {k:self.activation(v) for k,v in h.items()}

        with g.local_scope():
            g.ndata['h'] = h['node']
            hg = self.pool(g,h['node'])
            return self.classify(hg)

##################################################################################






with open('data/het_Dataset.pkl','rb') as f:
    g = pkl.load(f)
with open('data/het_train_ind.pkl','rb') as f:
    train_ind = pkl.load(f)
with open('data/het_valid_ind.pkl','rb') as f:
    valid_ind = pkl.load(f)
with open('data/het_test_ind.pkl','rb') as f:
    test_ind = pkl.load(f)
with open('data/het_labels.pkl','rb') as f:
    labels = pkl.load(f)

labels=torch.tensor(labels)
features = g.nodes['node'].data['feat']


def postProcScores(tp,fp,fn):
    nbClasses= len(tp)
    prec=[0 for i in range(nbClasses)]
    rec =[0 for i in range(nbClasses)]
    f1  =[0 for i in range(nbClasses)]

    for k in range(nbClasses):
        if tp[k]==0:
            continue
        prec[k] = tp[k]/(tp[k]+fp[k])
        rec[k] = tp[k]/(tp[k]+fn[k])
        f1[k] =2*(prec[k]*rec[k])/(prec[k]+rec[k])
    return prec,rec,f1

def evaluate(model, features, labels, mask):
    model.eval()
    with torch.no_grad():
        logits = model({'node':features})
        logits = logits[mask]
        labels = labels[mask]
        _, indices = torch.max(logits, dim=1)
        correct = torch.sum(indices == labels)
        return correct.item() * 1.0 / len(labels)

def scores(model, features, labels, mask,classes):
    model.eval()
    tp = [0 for i in range(classes)]
    fp = [0 for i in range(classes)]
    fn = [0 for i in range(classes)]
    with torch.no_grad():

        logits = model({'node':features})
        logits = logits[mask]
        labels = labels[mask]
        _, indices = torch.max(logits, dim=1)
        for one_cl in range(classes):
            tp[one_cl] = sum([1 for i in range(len(labels)) if (indices[i]==one_cl) and (labels[i]==one_cl)])
            fp[one_cl] = sum([1 for i in range(len(labels)) if (indices[i]==one_cl) and (labels[i]!=one_cl)])
            fn[one_cl] = sum([1 for i in range(len(labels)) if (indices[i]!=one_cl) and (labels[i]==one_cl)])

        return tp,fp,fn

loss_fcn = torch.nn.CrossEntropyLoss()

# use optimizer






for ite_glob in range(nbItes):
    tinit= time.time()
    model = Classifier(g,g.etypes,nbFeats,nbHidden,nbClasses,nbLayers,'sum',F.relu,mpool,drop)

    optimizer = torch.optim.Adam(model.parameters(),
                                 lr=lr,
                                 weight_decay=weight_decay)
    train_mask = torch.tensor([True if i in train_ind[ite_glob] else False for i in range(len(labels))])
    val_mask = torch.tensor([True if i in valid_ind[ite_glob] else False for i in range(len(labels))])
    test_mask = torch.tensor([True if i in test_ind[ite_glob] else False for i in range(len(labels))])

    nbTrain = sum([1 for i in range(len(labels)) if train_mask[i]])
    nbValid = sum([1 for i in range(len(labels)) if val_mask[i]])
    nbTest = sum([1 for i in range(len(labels)) if test_mask[i]])



    dur=[]
    liste_loss = []
    maxAcc=0
    for epoch in range(nbEpochs):
        model.train()
        if epoch >= 3:
            t0 = time.time()
            # forward
        logits = model({'node':features})
        # print(logits)
        loss = loss_fcn(logits[train_mask], labels[train_mask])
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        acc = evaluate(model, features, labels, val_mask)
        liste_loss.append(loss_fcn(logits[val_mask], labels[val_mask]).item())
        if epoch >= 3:
            dur.append(time.time() - t0)

            print("Epoch {:05d} | Time(s) {:.4f} | Loss {:.4f} | Accuracy {:.4f} | ". format(epoch, np.mean(dur), loss.item(),acc))
        if epoch>earlyStop and liste_loss[-1] > np.mean(liste_loss[-earlyStop-1:-1]):
            print("Early Stopping")
            break

    tp,fp,fn = scores(model, features, labels, test_mask,nbClasses)

    prec,rec,f1 = postProcScores(tp,fp,fn)
    print('prec : ',prec)
    print('rec : ',rec)
    print('f1 : ',f1)
    dicoScores['tp'].append(tp)
    dicoScores['fp'].append(fp)
    dicoScores['fn'].append(fn)
    dicoScores['stop'].append(epoch)
    dicoScores['time'].append(time.time()-tinit)

print(dicoScores)
with open(file2Save,'wb') as f:
    pkl.dump(dicoScores,f)
