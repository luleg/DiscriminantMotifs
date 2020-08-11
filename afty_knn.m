
addpath('Stock');
addpath('clusterix');
clear variables;clc
close all
load my_embs
addpath('clusterix')
nb_reseaux = size(V,2);

D = sqdist(V,V);
dmax = max(max(D));

sigma = 1/2*dmax/(nb_reseaux)^(1/9);
D = -1/sigma*D;
D = exp(D);
B = D;
D = D-diag(diag(D));

k = 20;
C = zeros(size(D));
for i=1:nb_reseaux
    ld = D(i,:);
    [ld,ind] = sort(ld,'descend');
    C(i,ind(1:k)) = ld(1:k);
end
[~,~,delta] = find(C); delta = range(delta);
figure(41),clf
spyc(C,5);
ttl = title('Closest Neighb. Sparsification');
ttl.FontSize= 20;


A = (C>0).*(9/delta*C+(10*delta-9)/delta);

[p] = clusterix(A,4,3,200,1);

%% Display Confusion matrix

clust = unique(p);

N = zeros(5,length(clust));
for i=1:length(classes)-1
    for j=1:length(clust)
        N(i,j) = sum(p(classes(i):classes(i+1)-1) == clust(j));
    end
end

Conf_knn = ["Food Webs";"Electronic Circuits";"Discourse Structures";"Social Networks";"LFR Networks"];
T = table(Conf_knn);
for j=1:length(clust)
    eval(['Cluster_',int2str(j),'=N(:,j);'])
    Tmp = eval(['table(Cluster_',int2str(j),')']);
    T = [T Tmp];
end
disp(T);

