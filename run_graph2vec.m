% script to create the k-nn classifier based on the discriminant motifs of
% the networks.
addpath('Stock')
clear variables;clc;
close all

addpath('..')

load RandPermDatasets
nb_ites = size(lind_elec,1);
% nb_ites = 10;

% some global variables

global disp_fig
disp_fig = 0; % do you want to display the Figures ?


% Parameters of the graph2vec algo 
% We have chosen those that provide the best scores
global t_emb
global t_wl
t_emb = 64;
t_wl = 1;



nb_basis = 40; % number of networks from each field used to generate the
% PCA basis
% (/!\ the smallest database is the electronic circuits : 52 networks)
% nb_motifs =9; % number of most discriminant motifs we want to take into
% account

tens_Conf = zeros(4,4,nb_ites); % confusion matrix:
% tens_Conf(i,j,k) : number of networks from field i that has been
% classified in field j by the classifier at the kth iteration


% Loading of the networks
load('names_fw')
nb_fw = length(lreseaux); % number of foodwebs
lreseaux_fw = lreseaux; % struc that contains the name of the networks

load('names_elec')
nb_elec = length(lreseaux);
lreseaux_elec = lreseaux;

load('names_stac')
nb_stac = length(lreseaux);
lreseaux_stac = lreseaux;

load('names_soc')
nb_soc = length(lreseaux);
lreseaux_soc = lreseaux;

nb_reseaux = nb_fw+nb_elec+nb_stac+nb_soc;

lreseaux = cell(1,nb_fw);
V_fw = zeros(t_emb,nb_fw);
for i=1:nb_fw
    load(['Graph2Vec/EmbMat/Wl',int2str(t_wl),'_Emb',int2str(t_emb),'/',lreseaux_fw{i}])
    V_fw(:,i) = emb; 
end

V_elec = zeros(t_emb,nb_elec);
for i=1:nb_elec
   load(['Graph2Vec/EmbMat/Wl',int2str(t_wl),'_Emb',int2str(t_emb),'/',lreseaux_elec{i}])
   V_elec(:,i) = emb; 
end

V_stac = zeros(t_emb,nb_stac);
for i=1:nb_stac
   load(['Graph2Vec/EmbMat/Wl',int2str(t_wl),'_Emb',int2str(t_emb),'/',lreseaux_stac{i}])
   V_stac(:,i) = emb; 
end

V_soc = zeros(t_emb,nb_soc);
for i=1:nb_soc
   load(['Graph2Vec/EmbMat/Wl',int2str(t_wl),'_Emb',int2str(t_emb),'/',lreseaux_soc{i}])
   V_soc(:,i) = emb; 
end


fw  = 1:nb_basis;
elec = fw+nb_basis;
stac = fw+2*nb_basis;
soc = fw+3*nb_basis;

% For the test stage
nbt_fw = nb_fw-nb_basis;
nbt_elec = nb_elec-nb_basis;
nbt_stac = nb_stac-nb_basis;
nbt_soc = nb_soc-nb_basis;
nbt_reseaux = nbt_fw+nbt_elec+nbt_stac+nbt_soc;
fwt  = 1:nbt_fw;
elect= 1:nbt_elec;elect= elect+nbt_fw;
stact= 1:nbt_stac;stact= stact+nbt_elec+nbt_fw;
soct = 1:nbt_soc; soct = soct +nbt_stac+nbt_elec+nbt_fw;


load RandPermDatasets

for iter = 1:nb_ites
    if mod(iter,50) == 0
        disp(['iteration number: ',int2str(iter),' over ',int2str(nb_ites),'.']);
    end
    
    
    ind_fw = lind_fw(iter,:);
    ind_elec = lind_elec(iter,:);
    ind_stac = lind_stac(iter,:);
    ind_soc = lind_soc(iter,:);
    
    C = [V_fw(:,ind_fw(1:nb_basis)),V_elec(:,ind_elec(1:nb_basis)),V_stac(:,ind_stac(1:nb_basis)),V_soc(:,ind_soc(1:nb_basis))];
 
    
    C = C';
    c_fw = 1/nb_basis*ones(1,nb_basis)*C(fw,:);
    c_elec = 1/nb_basis*ones(1,nb_basis)*C(elec,:);
    c_stac = 1/nb_basis*ones(1,nb_basis)*C(stac,:);
    c_soc = 1/nb_basis*ones(1,nb_basis)*C(soc,:);
    %% classif
    
    % Generation of the test basis: taking all the networks that have not
    % been used to generate the PCA basis, we will classify them with the
    % distance2mean classifier
    
    Ct = [V_fw(:,ind_fw(nb_basis+1:end)),V_elec(:,ind_elec(nb_basis+1:end)),V_stac(:,ind_stac(nb_basis+1:end)),V_soc(:,ind_soc(nb_basis+1:end))];

    Ct = Ct';
    dist2mean
    
    tens_Conf(:,:,iter) = Conf2mean;
end

Labels = ["Food Webs";"Electronic Circuits";"Discourse Structures";"Social Networks"];
[tab_prec,tab_rec,tab_f1,Conf] = postproc_scores(tens_Conf);

Precision = Labels;
Minimum = tab_prec(1,:)';
Maximum = tab_prec(2,:)';
Mean = tab_prec(3,:)';
Std = tab_prec(4,:)';
tab_prec = table(Precision,Minimum,Maximum,Mean,Std);

Recall = Labels;
Minimum = tab_rec(1,:)';
Maximum = tab_rec(2,:)';
Mean = tab_rec(3,:)';
Std = tab_rec(4,:)';
tab_rec = table(Recall,Minimum,Maximum,Mean,Std);

F1Score = Labels;
Minimum = tab_f1(1,:)';
Maximum = tab_f1(2,:)';
Mean = tab_f1(3,:)';
Std = tab_f1(4,:)';
tab_f1= table(F1Score,Minimum,Maximum,Mean,Std);

disp(tab_prec)
disp(tab_rec)
disp(tab_f1)