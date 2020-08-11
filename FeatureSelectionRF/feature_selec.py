import math
import numpy as np
import sklearn.metrics as sk
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from mat4py import loadmat

def recup_features(name_file):
    vector = loadmat(folder_to_matFiles+'MatNetworks/'+name_file+'.mat')
    vector = vector['Pbm']
    mot3 = vector['motif3']
    mot3.sort(key = lambda a:a[0])
    mot3 = [elt[1] for elt in mot3]
    mot4 = vector['motif4']
    mot4.sort(key = lambda a:a[0])
    mot4 = [elt[1] for elt in mot4]
    vector = mot3+mot4
    return vector

def recup_features_preproc(name_file):
    vector = loadmat(folder_to_matFiles+'MatNetworks/'+name_file+'.mat')
    vector = vector['Pbm']
    nb_nodes = vector['nb_nodes']
    mot3 = vector['motif3']
    mot3.sort(key = lambda a:a[0])
    mot3 = np.array([elt[1] for elt in mot3])
    mot4 = vector['motif4']
    mot4.sort(key = lambda a:a[0])
    mot4 = np.array([elt[1] for elt in mot4])
    mot3 = 13/(nb_nodes*(nb_nodes-1)*(nb_nodes-2))*mot3
    mot4 = 199/(nb_nodes*(nb_nodes-1)*(nb_nodes-2)*(nb_nodes))*mot4
    vector = np.concatenate((mot3,mot4))
    vector = vector/np.linalg.norm(vector)
    return vector


folder_to_matFiles = '../Stock/'

vector = loadmat(folder_to_matFiles+'MatNetworks/elec_b03_0.mat')

vector = vector['Pbm']
nb_nodes = vector['nb_nodes']
mot3 = vector['motif3']
mot3.sort(key = lambda a:a[0])
id3 = [elt[0] for elt in mot3]
mot4 = vector['motif4']
mot4.sort(key = lambda a:a[0])
id4 = [elt[0] for elt in mot4]


names_elec = loadmat(folder_to_matFiles+'names_elec.mat')
names_elec = names_elec['lreseaux']
names_fw = loadmat(folder_to_matFiles+'names_fw.mat')
names_fw = names_fw['lreseaux']
names_soc = loadmat(folder_to_matFiles+'names_soc.mat')
names_soc = names_soc['lreseaux']
names_stac = loadmat(folder_to_matFiles+'names_stac.mat')
names_stac = names_stac['lreseaux']

randperm = loadmat(folder_to_matFiles+'RandPermDatasets.mat')
randperm_elec = randperm['lind_elec']
randperm_fw = randperm['lind_fw']
randperm_soc = randperm['lind_soc']
randperm_stac = randperm['lind_stac']

lsel_features = []
lprec = []
lrec = []
lf1 = []

nb_test = 40
nb_ites = 500
for k in range(nb_ites):
    X_train = []
    y_train = []
    X_test = []
    y_test = []
    for ind in randperm_elec[k][:nb_test]:
        # load and process the feature vector
        crt_file = names_elec[ind-1]
        vector = recup_features_preproc(crt_file)
        X_train.append(vector)
        y_train.append(1)

    for ind in randperm_fw[k][:nb_test]:
        crt_file = names_fw[ind-1]
        vector = recup_features_preproc(crt_file)
        X_train.append(vector)
        y_train.append(0)

    for ind in randperm_soc[k][:nb_test]:
        crt_file = names_soc[ind-1]
        vector = recup_features_preproc(crt_file)
        X_train.append(vector)
        y_train.append(3)

    for ind in randperm_stac[k][:nb_test]:
        crt_file = names_stac[ind-1]
        vector = recup_features_preproc(crt_file)
        X_train.append(vector)
        y_train.append(2)

    model = RandomForestClassifier(n_estimators = 100,random_state=0).fit(X_train,y_train)
    lsel_features.append(model.feature_importances_)

    print("end of iteration : ",k)


str_sel_features = ""
for sel_features in lsel_features:
    for elt in sel_features:
        str_sel_features+=str(elt)+','
    str_sel_features+='\n'
with open('selected_features.txt','w') as f:
    f.write(str_sel_features)
