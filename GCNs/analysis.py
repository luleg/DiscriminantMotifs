import os
import pickle as pkl
import argparse as arg
import numpy as np


parser = arg.ArgumentParser(description="Analyse the quality of the (R)GCN classifier, whose scores are stored in a pickle dict.\n-file: path o the pickl file containing the dict.")
parser.add_argument('-file',required=True)

args = parser.parse_args()

with open(args.file,'rb') as f:
    scores = pkl.load(f)


tps = scores['tp']
fps = scores['fp']
fns = scores['fn']
nbIte = len(tps)

prec = [[0,0,0,0] for k in range(len(tps))]
rec = [[0,0,0,0] for k in range(len(tps))]
f1 = [[0,0,0,0] for k in range(len(tps))]
prec2 = [[],[],[],[]]
rec2 = [[],[],[],[]]
f12 = [[],[],[],[]]
ite_nul=False
lite = []
times = 0
epochs = 0
for i in range(nbIte):
    tp = tps[i]
    fp = fps[i]
    fn = fns[i]
    for k in range(4):
        if tp[k] == 0:
            ite_nul=True
            prec2[k].append(0)
            rec2[k].append(0)
            f12[k].append(0)
            continue
        prec[i][k] =tp[k]/(tp[k]+fp[k])
        rec[i][k] =tp[k]/(tp[k]+fn[k])
        f1[i][k] = 2*(prec[i][k]*rec[i][k])/(prec[i][k]+rec[i][k])
        prec2[k].append(prec[i][k])
        rec2[k].append(rec[i][k])
        f12[k].append(f1[i][k])
    if ite_nul:
        print('Empty Class: --> ',scores['stop'][i],'--',scores['time'][i])
        lite.append(ite_nul)
        ite_nul=False
    else:
        # print('OK:-->',scores['stop'][i],'--',scores['time'][i])
        lite.append(ite_nul)
        times+=scores['time'][i]
        epochs+=scores['stop'][i]


mprec = [np.mean([l for it,l in enumerate(ll) if not lite[it]]) for ll in prec2]
mrec = [np.mean([l for it,l in enumerate(ll) if not lite[it]]) for ll in rec2]
mf1 = [np.mean([l for it,l in enumerate(ll) if not lite[it]]) for ll in f12]

sprec = [np.std(ll) for ll in prec2]
srec = [np.std(ll) for ll in rec2]
sf1 = [np.std(ll) for ll in f12]

print('Without Empty Classes -------------------------------------------')
print('Precision : ')
print('\t'.join(['{:.4f}'.format(k) for k in mprec]))
print('\t'.join(['{:.4f}'.format(k) for k in sprec]),'\n')

print('Recall : ')
print('\t'.join(['{:.4f}'.format(k) for k in mrec]))
print('\t'.join(['{:.4f}'.format(k) for k in srec]),'\n')

print('F1 : ')
print('\t'.join(['{:.4f}'.format(k) for k in mf1]))
print('\t'.join(['{:.4f}'.format(k) for k in sf1]),'\n')

print('Nb non considered : ',len([i for i in lite if i]))

print('With Empty Classes-------------------------------------------')

mprec = [np.mean([l for it,l in enumerate(ll)]) for ll in prec2]
mrec = [np.mean([l for it,l in enumerate(ll)]) for ll in rec2]
mf1 = [np.mean([l for it,l in enumerate(ll)]) for ll in f12]

sprec = [np.std(ll) for ll in prec2]
srec = [np.std(ll) for ll in rec2]
sf1 = [np.std(ll) for ll in f12]

print('Precision : ')
print('\t'.join(['{:.4f}'.format(k) for k in mprec]))
print('\t'.join(['{:.4f}'.format(k) for k in sprec]),'\n')

print('Recall : ')
print('\t'.join(['{:.4f}'.format(k) for k in mrec]))
print('\t'.join(['{:.4f}'.format(k) for k in srec]),'\n')

print('F1 : ')
print('\t'.join(['{:.4f}'.format(k) for k in mf1]))
print('\t'.join(['{:.4f}'.format(k) for k in sf1]),'\n')

print("-----------------------------------------------------------------------")
print("Time : ",times,"s")
print("Number of epochs: ",epochs," for ",nbIte," tests.")
# [0.6207362640799174, 0.8137993601228896, 0.9874727668845316, 0.9715927973984644]
# [0.5555555555555557, 0.8962962962962963, 0.8238482384823849, 1.0]
# [0.5658808652865883, 0.8511351664126098, 0.89794175939814, 0.9855552541332577]
