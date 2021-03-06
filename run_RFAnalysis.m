addpath('Stock')
addpath('FeatureSelectionRF')
clear variables;clc;
close all

features_RF = csvread('selected_features.txt');
features_RF(:,end)= [];
disp_gamma_errorbar(features_RF,51,30,true,0,-1,'RF-Score')


% script to create the k-nn classifier based on the discriminant motifs of
% the networks.
global disp_fig
disp_fig = false; % do you want to display the Figures ?

load RandPermDatasets
nb_ites = size(lind_elec,1);
% nb_ites = 10;


nb_basis = 40; % number of networks from each field used to generate the
% PCA basis
% (/!\ the smallest database is the electronic circuits : 52 networks)
% nb_motifs =9; % number of most discriminant motifs we want to take into
% account
global nb_motifs
nb_motifs = 9;

global keep_motifs
global t_vect
t_vect = 13+199;
% t_vect = nb_motifs;
% t_vect = 199+13; % size of the vectors used to represent networks :
% 199 4-node motifs, 13 3-node motifs

% Parametrisation of the 2-step normalisation  :
global type_norm
type_norm = 0; % first step
% 0 : 13*v3/(n-2)*(n-1)*n
% 1 : 13*6*v3/(n-2)*(n-1)*n
% 2 : 13*v3/norm(v3), 199*v4/norm(v4)
global norm_2
norm_2 = 0; % second step (normalisation in norm-2)
% 1 : x/norm(x)
% 2 : x(motifk)/norm([xi(motifk), i = 1 :nb_graphs])
% else :  no normalisation


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
for i=1:nb_fw
    lreseaux{i} = lreseaux_fw{i};
end
V_fw = generate_basis(lreseaux);

lreseaux = cell(1,nb_elec);
for i=1:nb_elec
    lreseaux{i} = lreseaux_elec{i};
end
V_elec = generate_basis(lreseaux);

lreseaux = cell(1,nb_stac);
for i=1:nb_stac
    lreseaux{i} = lreseaux_stac{i};
end
V_stac = generate_basis(lreseaux);

lreseaux = cell(1,nb_soc);
for i=1:nb_soc
    lreseaux{i} = lreseaux_soc{i};
end
V_soc = generate_basis(lreseaux);


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
    Gini_crt = features_RF(iter,:);
    [~,sort_motif] = sort(Gini_crt,'descend');
    keep_motifs = sort_motif(1:nb_motifs);
    
    
    ind_fw = lind_fw(iter,:);
    ind_elec = lind_elec(iter,:);
    ind_stac = lind_stac(iter,:);
    ind_soc = lind_soc(iter,:);
    
    C = [V_fw(:,ind_fw(1:nb_basis)),V_elec(:,ind_elec(1:nb_basis)),V_stac(:,ind_stac(1:nb_basis)),V_soc(:,ind_soc(1:nb_basis))];
    
    % postproc :
    C = C(keep_motifs,:);
    N = ones(1,length(keep_motifs))*(C.*C);N = diag(1./sqrt(N));
    C = C*N;
    
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
    Ct = Ct(keep_motifs,:);
    N = ones(1,length(keep_motifs))*(Ct.*Ct);N = diag(1./sqrt(N));
    Ct = Ct*N;
    
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