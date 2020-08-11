
addpath('Stock');
addpath('clusterix');
clear variables;clc;
close all
load my_embs
addpath('clusterix')

nb_reseaux = size(V,2);

%% OK for sparse-threshold. If does not work change tol in EdgeRefinement (line 27 -> should be 3)
D = sqdist(V,V);
dmax = max(max(D));

sigma = 1/2*dmax/(nb_reseaux)^(1/9);
D = -1/sigma*D;
D = exp(D);
B = D;
D = D-diag(diag(D));
figure(12),clf,
spyc(D,5)

ttl = title('Affinity Matrix');
ttl.FontSize = 20;


A = D.*(D>5e-1);
figure(51),clf,
spyc(A,5)
ttl = title('Threshold Sparsification');
ttl.FontSize = 20;




A = A-diag(diag(A));
[~,~,delta] = find(A); delta = range(delta);
A = (A>0).*(9/delta*A+(10*delta-9)/delta);
A = 1/2*(A+A');


[p] = clusterix(A,3,5,11,0);


%% Display Confusion Matrix
clust = unique(p);

N = zeros(5,length(clust));

for i=1:length(classes)-1
    for j=1:length(clust)
        N(i,j) = sum(p(classes(i):classes(i+1)-1) == clust(j));
    end
end

Conf_Threshold = ["Food Webs";"Electronic Circuits";"Discourse Structures";"Social Networks";"LFR Networks"];
T = table(Conf_Threshold);
for j=1:length(clust)
    eval(['Cluster_',int2str(j),'=N(:,j);'])
    Tmp = eval(['table(Cluster_',int2str(j),')']);
    T = [T Tmp];    
end
disp(T);
