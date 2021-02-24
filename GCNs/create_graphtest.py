import dgl
import torch
import networkx as nx
from collections import Counter
import pickle as pkl


def LoadDatasetSimple(path2lists):
    # path2lists = '../../Documents/GitHub/DiscriminantMotifs/Stock/'
    elec = 'elec'
    fw = 'fw'
    soc= 'soc'
    stac='stac'
    classes_list = [elec,fw,soc,stac]
    dico_lab = {classes_list[i]:i for i in range(len(classes_list))}

    path2nwks = 'Processing/Networks/'


    ##
    dico_name_nwks = {one_class:[] for one_class in classes_list}
    # print(dico_nwks)
    for one_class in classes_list:
            with open(path2lists+'liste_'+one_class+'.txt','r') as f:
                for line in f:
                    line=line.split()
                    dico_name_nwks[one_class].append(line[0])


    dico_nwks = {one_class:[] for one_class in classes_list}
    dico_name = dict()
    liste_4bgl = []
    labels=[]
    num_nwk=-1
    for one_class,liste_nwks in dico_name_nwks.items():
        print('Create homogeneous networks of type: ',one_class)

        for one_nwk in liste_nwks:
            num_nwk+=1
            # print(path2lists+path2nwks+one_nwk)
            dico_name[one_nwk]=num_nwk
            G = nx.read_edgelist(path2lists+path2nwks+one_nwk+'.nwk',nodetype=int,create_using=nx.DiGraph)
            nbNodes= G.order()
            G.add_edges_from([(i,i) for i in range(nbNodes)])
            g = dgl.from_networkx(G)
            features=[[float(G.in_degree[i]),float(G.out_degree[i])] for i in range(nbNodes)]
            g.ndata['h'] = torch.tensor(features)
            liste_4bgl.append(g)

            labels.append(dico_lab[one_class])
    bg = dgl.batch(liste_4bgl)

    train_msk= [[] for i in range(500)]
    valid_msk= [[] for i in range(500)]
    test_msk = [[] for i in range(500)]
    nbTrain=35
    nbValid=5

    for one_class in classes_list:
        with open(path2lists+'perms_'+one_class+'.txt','r') as f:
            num_ite = 0
            for line in f:
                nbTrainTmp=0
                nbValidTmp=0
                line=line.split()
                if len(line) == 0:
                    continue
                line = line[0].split(',')
                for nwk in line:
                    if (nbTrainTmp<nbTrain):
                        train_msk[num_ite].append(dico_name[nwk])
                        nbTrainTmp+=1
                    elif (nbValidTmp<nbValid):
                        valid_msk[num_ite].append(dico_name[nwk])
                        nbValidTmp+=1
                num_ite +=1

    for num_ite in range(500):
        for i in range(len(labels)):
            if (i not in train_msk[num_ite]) and (i not in valid_msk[num_ite]):
                test_msk[num_ite].append(i)


    return bg,labels,train_msk,valid_msk,test_msk



def LoadDatasetHetero(path2lists):
    # path2lists = '../../Documents/GitHub/DiscriminantMotifs/Stock/'
    elec = 'elec'
    fw = 'fw'
    soc= 'soc'
    stac='stac'
    classes_list = [elec,fw,soc,stac]
    dico_lab = {classes_list[i]:i for i in range(len(classes_list))}




    path2nwks = 'Processing/Networks/'


    ##
    dico_name_nwks = {one_class:[] for one_class in classes_list}
    # print(dico_nwks)
    for one_class in classes_list:
            with open(path2lists+'liste_'+one_class+'.txt','r') as f:
                for line in f:
                    line=line.split()
                    dico_name_nwks[one_class].append(line[0])


    dico_nwks = {one_class:[] for one_class in classes_list}
    dico_name = dict()
    liste_4bgl = []
    labels=[]
    num_nwk=-1
    for one_class,liste_nwks in dico_name_nwks.items():
        print('Create heterogeneous networks of type: ',one_class)

        for one_nwk in liste_nwks:
            num_nwk+=1
            # print(path2lists+path2nwks+one_nwk)
            dico_name[one_nwk]=num_nwk
            with open(path2lists+path2nwks+one_nwk+'.nwk','r') as f:
                contenu = f.read()
            contenu = contenu.split('!m:')[-1]
            contenu = list(filter(lambda a:a!='',contenu.split('\n')[1:]))
            contenu = [(int(line.split()[0]),int(line.split()[1])) for line in contenu]
            lsrcs = [u[0] for u in contenu]
            ltgts = [u[1] for u in contenu]

            nb_node = max(lsrcs+ltgts) + 1
            cdin = Counter(lsrcs)
            cdout = Counter(ltgts)

            lsrcs += [i for i in range(nb_node)]
            ltgts += [i for i in range(nb_node)]


            g = dgl.heterograph({
                ('node','points_to','node'):(lsrcs,ltgts),
                ('node','pointed_by','node'):(ltgts,lsrcs)
            })
            features = [[float(cdin[i]),float(cdout[i])] for i in range(nb_node)]
            g.nodes['node'].data['feat'] = torch.tensor(features)

            liste_4bgl.append(g)

            labels.append(dico_lab[one_class])
    bg = dgl.batch(liste_4bgl)

    train_msk= [[] for i in range(500)]
    valid_msk= [[] for i in range(500)]
    test_msk = [[] for i in range(500)]
    nbTrain=35
    nbValid=5

    for one_class in classes_list:
        with open(path2lists+'perms_'+one_class+'.txt','r') as f:
            num_ite = 0
            for line in f:
                nbTrainTmp=0
                nbValidTmp=0
                line=line.split()
                if len(line) == 0:
                    continue
                line = line[0].split(',')
                for nwk in line:
                    if (nbTrainTmp<nbTrain):
                        train_msk[num_ite].append(dico_name[nwk])
                        nbTrainTmp+=1
                    elif (nbValidTmp<nbValid):
                        valid_msk[num_ite].append(dico_name[nwk])
                        nbValidTmp+=1
                num_ite +=1

    for num_ite in range(500):
        for i in range(len(labels)):
            if (i not in train_msk[num_ite]) and (i not in valid_msk[num_ite]):
                test_msk[num_ite].append(i)


    return bg,labels,train_msk,valid_msk,test_msk





if __name__ == '__main__':
    path2lists = '../Stock/'
    bg,labels,train_msk,valid_msk,test_msk = LoadDatasetSimple(path2lists)

    with open('data/Dataset.pkl','wb') as f:
        pkl.dump(bg,f)
    with open('data/train_ind.pkl','wb') as f:
        pkl.dump(train_msk,f)
    with open('data/valid_ind.pkl','wb') as f:
        pkl.dump(valid_msk,f)
    with open('data/test_ind.pkl','wb') as f:
        pkl.dump(test_msk,f)
    with open('data/labels.pkl','wb') as f:
        pkl.dump(labels,f)


    bg,labels,train_msk,valid_msk,test_msk = LoadDatasetHetero(path2lists)

    with open('data/het_Dataset.pkl','wb') as f:
        pkl.dump(bg,f)
    with open('data/het_train_ind.pkl','wb') as f:
        pkl.dump(train_msk,f)
    with open('data/het_valid_ind.pkl','wb') as f:
        pkl.dump(valid_msk,f)
    with open('data/het_test_ind.pkl','wb') as f:
        pkl.dump(test_msk,f)
    with open('data/het_labels.pkl','wb') as f:
        pkl.dump(labels,f)
