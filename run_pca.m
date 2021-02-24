% script to create the distance-to-mean classifier based on the PCA of the
% networks in the training set, decomposed onto 3-and/or 4-and/or 5-node 
% graphlets.
addpath('Stock')
clear variables;clc;
close all

% some global variables

% k-node graphlets to use:
% motif3 = 1; 
% motif4 = 1; 
% motif5 = 0; 
% means 3-and-4-node graphlets are used within the analysis, but not 5-node
% graphlets 
motif3 = 1;
motif4 = 1;
motif5 = 0;

if motif5
    load RandPermDatasets5
else
    load RandPermDatasets
end
% nb_ites = size(lind_elec,1); 
nb_ites = 500;

global disp_fig
disp_fig = 0; % do you want to display the Figures ?

global p
p =-1; % percentage of trace to keep in the PCA (if <0 one keeps a number
% of axes instead of a percentage of trace)
global nb_ax
nb_ax = 10;

nb_basis = 40; % number of networks from each field used to generate the
% PCA basis
% (/!\ the smallest database is the electronic circuits : 52 networks)

% Parametrisation of the 2-step normalisation  :
global type_norm
type_norm = 0; % first step
% 0 : alpha_k = m_k/k!/(k among n)
% 1 : alpha_k = m_k/(k among n)
% 2 : alpha_k = 1/(k among n)
% 3 : alpha_k = log(mk)/k!/(k among n)
% 4 : alpha_k = log(mk/k!/(k among n))
% 5 : alpha_k = 1/||Phi_k||
% 6 : alpha_k = sqrt(m_k)/||Phi_k||
% 7 : alpha_k = 1
global norm_2
norm_2 = 1; % second step (normalisation in norm-2)
% 1 : x/norm(x)
% 2 : x(motifk)/norm([xi(motifk), i = 1 :nb_graphs])
% else :  no normalisation



tens_per = zeros(1,nb_ites); % number of ppal axes to obtain the target %
% of trace (or % of trce given a fie number of ppal axes)
tens_Conf = zeros(4,4,nb_ites);
% tens_Conf(i,j,k) : number of networks from field i that has been
% classified in field j by the classifier at the kth iteration


t_genere_basis = cputime;
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

if motif5
    load('names_soc5')
else
    load('names_soc')
end
nb_soc = length(lreseaux);
lreseaux_soc = lreseaux;

nb_reseaux = nb_fw+nb_elec+nb_stac+nb_soc;

lreseaux = cell(1,nb_fw);
for i=1:nb_fw
    lreseaux{i} = lreseaux_fw{i};
end
V_fw = generate_basis(lreseaux,motif3,motif4,motif5);

lreseaux = cell(1,nb_elec);
for i=1:nb_elec
    lreseaux{i} = lreseaux_elec{i};
end
V_elec = generate_basis(lreseaux,motif3,motif4,motif5);

lreseaux = cell(1,nb_stac);
for i=1:nb_stac
    lreseaux{i} = lreseaux_stac{i};
end
V_stac = generate_basis(lreseaux,motif3,motif4,motif5);

lreseaux = cell(1,nb_soc);
for i=1:nb_soc
    lreseaux{i} = lreseaux_soc{i};
end
V_soc = generate_basis(lreseaux,motif3,motif4,motif5);

t_genere_basis = cputime-t_genere_basis;

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

t_train = 0;
t_test = 0;

% For each iteration
for iter = 1:nb_ites
    if mod(iter,10) == 0
        disp(['iteration number: ',int2str(iter),' over ',int2str(nb_ites),'.']);
    end
    ind_fw = lind_fw(iter,:);
    ind_elec = lind_elec(iter,:);
    ind_stac = lind_stac(iter,:);
    ind_soc = lind_soc(iter,:);
    temps = cputime;
    %
    %% PCA
    % Generation of the basis
       
    V = [V_fw(:,ind_fw(1:nb_basis)),V_elec(:,ind_elec(1:nb_basis)),V_stac(:,ind_stac(1:nb_basis)),V_soc(:,ind_soc(1:nb_basis))];
    
    
    % calling of the acp script, that manages the PCA process
    pca
    tens_per(iter)=prct_trace; % Knowing percentage of trace that is kept
    
    % mean individuals for each benchmark :
    c_fw = 1/nb_basis*ones(1,nb_basis)*C(fw,:);
    c_elec = 1/nb_basis*ones(1,nb_basis)*C(elec,:);
    c_stac = 1/nb_basis*ones(1,nb_basis)*C(stac,:);
    c_soc = 1/nb_basis*ones(1,nb_basis)*C(soc,:);
    t_train = t_train + cputime-temps;
    %% classif
    
    % Generation of the test basis: taking all the networks that have not
    % been used to generate the PCA basis, we will classify them with the
    % distance2mean classifier
    temps= cputime;
    
    Vt =[V_fw(:,ind_fw(nb_basis+1:end)),V_elec(:,ind_elec(nb_basis+1:end)),V_stac(:,ind_stac(nb_basis+1:end)),V_soc(:,ind_soc(nb_basis+1:end))];
    % change of frame for the test dataset according to the PCA basis
    Ct = (Vt-x_b*ones(1,nbt_reseaux))'*U(:,1:nb_ax);
    
    % call to classifier
    dist2mean
    t_test = t_test + cputime-temps;
    tens_Conf(:,:,iter) = Conf2mean;
    
end


Labels = ["Food Webs";"Electronic Circuits";"Discourse Structures";"Social Networks"];
[tab_prec,tab_rec,tab_f1,Conf,tab_acc] = postproc_scores(tens_Conf);

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

F1_Score = Labels;
Minimum = tab_f1(1,:)';
Maximum = tab_f1(2,:)';
Mean = tab_f1(3,:)';
Std = tab_f1(4,:)';
tab_f1= table(F1_Score,Minimum,Maximum,Mean,Std);


Accuracy = "Global Accuracy";
Minimum = tab_acc(1);
Maximum = tab_acc(2);
Mean = tab_acc(3);
StandardDev = tab_acc(4);

tab_acc= table(Accuracy,Minimum,Maximum,Mean,StandardDev);

disp(tab_prec)
disp(tab_rec)
disp(tab_f1)
disp(tab_acc)


