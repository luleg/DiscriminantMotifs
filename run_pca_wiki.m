% script to create the distance-to-mean classifier based on the PCA of the
% wikipedia networks.
addpath('Stock')
clear variables;clc;

close all
nb_ites = 10;

% some global variables
motif3 = 1;
motif4 = 0;
motif5 = 0;
global disp_fig
disp_fig = 0; % do you want to display the Figures ?

global p
p =-1; 
global nb_ax
nb_ax =4;

% Parametrisation of the 2-step normalisation  :
global type_norm
type_norm = 0; 
global norm_2
norm_2 = 1; 

t_genere_basis = cputime;
% Loading of the networks
load('names_wiki_qualite')
nb_quali = length(lreseaux); % number of quality articles
lreseaux_quali = lreseaux; % struc that contains the name of the networks

load('names_wiki_non-neutre')
nb_nn = length(lreseaux);
lreseaux_nn = lreseaux;

load('names_wiki_guerre-edition')
nb_ge = length(lreseaux);
lreseaux_ge = lreseaux;

nb_reseaux = nb_quali+nb_nn+nb_ge;

lreseaux = cell(1,nb_quali);
for i=1:nb_quali
    lreseaux{i} = lreseaux_quali{i};
end
V_quali = generate_basis(lreseaux,motif3,motif4,motif5);

lreseaux = cell(1,nb_nn);
for i=1:nb_nn
    lreseaux{i} = lreseaux_nn{i};
end
V_nn= generate_basis(lreseaux,motif3,motif4,motif5);

lreseaux = cell(1,nb_ge);
for i=1:nb_ge
    lreseaux{i} = lreseaux_ge{i};
end
V_ge= generate_basis(lreseaux,motif3,motif4,motif5);

t_genere_basis = cputime-t_genere_basis;

V_pp = [V_nn,V_ge];
nbBasis = 20;

quali = 1:nbBasis;
confl = 1:nbBasis;confl = confl+nbBasis;

classesC = zeros(2,2,10);
for it= 1 :nb_ites
permsQ = randperm(33);
permsC = randperm(55);
V_quali = V_quali(:,permsQ);
V_pp= V_pp(:,permsC);

V = [V_quali(:,1:nbBasis),V_pp(:,1:nbBasis)];
Vt = [V_quali(:,nbBasis+1:end),V_pp(:,nbBasis+1:end)];
% V = [V_quali,V_pp];
pca

mu_q = 1/nbBasis*ones(1,nbBasis)*C(quali,:);
mu_c = 1/nbBasis*ones(1,nbBasis)*C(confl,:);

Ct = (Vt-x_b)'*U(:,1:nb_ax);

classes = fct_dist2mean([mu_q;mu_c]',Ct',[1,2]);
classesC(1,1,it) = sum(classes(1:13) == 1);
classesC(2,1,it) = sum(classes(1:13) == 2);
classesC(1,2,it) = sum(classes(14:end) == 1);
classesC(2,2,it) = sum(classes(14:end) == 2);
end

prec=zeros(2,10);
rec=zeros(2,10);
f1 = zeros(2,10);

for it=1:nb_ites
    prec(1,it) = classesC(1,1,it)/(classesC(1,1,it)+classesC(1,2,it));
    prec(2,it) = classesC(2,2,it)/(classesC(2,2,it)+classesC(2,1,it));
    rec(1,it) = classesC(1,1,it)/(classesC(1,1,it)+classesC(2,1,it));
    rec(2,it) = classesC(2,2,it)/(classesC(2,2,it)+classesC(1,2,it));
    f1(1,it) = 2*(prec(1,it)*rec(1,it))/(prec(1,it)+rec(1,it));
    f1(2,it) = 2*(prec(2,it)*rec(2,it))/(prec(2,it)+rec(2,it));
end

Labels=["Quality";"Contested"];
Mean = [mean(f1(1,:));mean(f1(2,:))];
stand_dev = [std(f1(1,:));std(f1(2,:))];
F1_score = table(Mean,stand_dev);

Mean = [mean(prec(1,:));mean(prec(2,:))];
stand_dev = [std(prec(1,:));std(prec(2,:))];
Precision = table(Mean,stand_dev);

Mean = [mean(rec(1,:));mean(rec(2,:))];
stand_dev = [std(rec(1,:));std(rec(2,:))];
Recall = table(Mean,stand_dev);

results = table(Labels,Precision,Recall,F1_score);
disp(results)