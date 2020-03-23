% script to create the k-nn classifier based on the discriminant motifs of
% the networks.
clear variables;clc;
close all

% some global variables

global nb_ites
nb_ites = 500; % number of iterations (number of classifiers based on 1 PCA
% basis instance)
global disp_fig
disp_fig = false; % do you want to display the Figures ?

nb_basis = 40; % number of networks from each field used to generate the
% PCA basis
% (/!\ the smallest database is the electronic circuits : 52 networks)
nb_motifs = 7; % number of most discriminant motifs we want to take into
% account
load('discriminant_motifs')
keep_motifs = sort_motif(1:nb_motifs,1);

% Parametrisation of the 2-step normalisation  :
global type_norm
type_norm = 0; % first step
% 0 : 13*v3/(n-2)*(n-1)*n
% 1 : 13*6*v3/(n-2)*(n-1)*n
% 2 : 13*v3/norm(v3), 199*v4/norm(v4)
global norm_2
norm_2 = true; % second step (normalisation in norm-2)

global nb_k
nb_k = 7;% for the k-nn classifier number of closest neighbours to look at

tens_infos = zeros(6,4,nb_ites); % scores to evaluate the classifier quality
% tens_info(1,k,i) : precision of kth field of networks at ite i
% tens_info(2,k,i) : recall of kth field of networks at ite i
% tens_info(3,k,i) : f1-score for kth field of networks at ite i
% tens_info(4,k,i) : true positive of kth field of networks at ite i
% tens_info(5,k,i) : fasle positive of kth field of networks at ite i
% tens_info(6,k,i) : false negative of kth field of networks at ite i

tens_Conf = zeros(4,4,nb_ites); % confusion matrix:
% tens_Conf(i,j,k) : number of networks from field i that has been
% classified in field j by the classifier at the kth iteration

global t_vect
t_vect = 199+13; % size of the vectors used to represent networks :
% 199 4-node motifs, 13 3-node motifs

% Loading of the networks
load('names_fw')
nb_fw = length(lreseaux); % number of foodwebs
lreseaux_fw = lreseaux; % struc that contains the name of the netwokrs

load('names_elec')
nb_elec = length(lreseaux);
lreseaux_elec = lreseaux;

load('names_stac')
nb_stac = length(lreseaux);
lreseaux_stac = lreseaux;

load('names_soc')
nb_soc = length(lreseaux);
lreseaux_soc = lreseaux;

% number of networks in each field (i for initial)
nb_fwi = nb_fw;
nb_eleci = nb_elec;
nb_staci = nb_stac;
nb_soci = nb_soc;

% Management of indices :
% For the learning stage
nb_fw = nb_basis;
nb_elec = nb_basis;
nb_stac = nb_basis;
nb_soc = nb_basis;
nb_reseaux = nb_fw+nb_elec+nb_stac+nb_soc;
fw  = 1:nb_fw;
elec= 1:nb_elec;elec= elec+nb_fw;
stac= 1:nb_stac;stac= stac+nb_elec+nb_fw;
soc = 1:nb_soc; soc = soc +nb_stac+nb_elec+nb_fw;
% For the test stage
nbt_fw = nb_fwi-nb_fw;
nbt_elec = nb_eleci-nb_elec;
nbt_stac = nb_staci-nb_stac;
nbt_soc = nb_soci-nb_soc;
nbt_reseaux = nbt_fw+nbt_elec+nbt_stac+nbt_soc;
fwt  = 1:nbt_fw;
elect= 1:nbt_elec;elect= elect+nbt_fw;
stact= 1:nbt_stac;stact= stact+nbt_elec+nbt_fw;
soct = 1:nbt_soc; soct = soct +nbt_stac+nbt_elec+nbt_fw;

% For each iteration
for iter = 1:nb_ites
    if mod(iter,50) == 0
        iter
    end
    ind_fw = randperm(nb_fwi);
    ind_elec = randperm(nb_eleci);
    ind_stac = randperm(nb_staci);
    ind_soc = randperm(nb_soci);
    
    lreseaux = {};
    
    for i=1:nb_fw
        lreseaux{i} = lreseaux_fw{ind_fw(i)};
    end
    for i=1:nb_elec
        lreseaux{i+nb_fw} = lreseaux_elec{ind_elec(i)};
    end
    for i=1:nb_stac
        lreseaux{i+nb_fw+nb_elec} = lreseaux_stac{ind_stac(i)};
    end
    for i=1:nb_soc
        lreseaux{i+nb_fw+nb_elec+nb_stac} = lreseaux_soc{ind_soc(i)};
    end
    
    V = generate_basis(lreseaux);
    
    
    %% knn
    
    % Generation of the test basis: taking all the networks that have not
    % been used to generate the PCA basis, we will classify them with a
    % k-nn classifier:
    lreseaux_test = {};
    
    for i=1:nbt_fw
        lreseaux_test{i} = lreseaux_fw{ind_fw(nb_fw+i)};
    end
    for i=1:nbt_elec
        lreseaux_test{i+nbt_fw} = lreseaux_elec{ind_elec(nb_elec+i)};
    end
    for i=1:nbt_stac
        lreseaux_test{i+nbt_fw+nbt_elec} = lreseaux_stac{ind_stac(nb_stac+i)};
    end
    for i=1:nbt_soc
        lreseaux_test{i+nbt_fw+nbt_elec+nbt_stac} = lreseaux_soc{ind_soc(nb_soc+i)};
    end
    
    Vt = generate_basis(lreseaux_test);
    
    C =V(keep_motifs,:)';
    Ct = Vt(keep_motifs,:)';
    
    knn
    
    V_all = [C;Ct];
    fw_all = [fw,nb_reseaux+fwt];nb_fw_all = length(fw_all);
    elec_all = [elec,nb_reseaux+elect];nb_elec_all = length(elec_all);
    stac_all = [stac,nb_reseaux+stact];nb_stac_all = length(stac_all);
    soc_all = [soc,nb_reseaux+soct];nb_soc_all = length(soc_all);
    
    
    tens_infos(:,:,iter,nb_motifs) = infos;
    tens_Conf(:,:,iter,nb_motifs) = Conf;
    
end

[tab_prec,tab_rec,tab_f1,Conf] = postproc_scores(tens_infos(:,:,:,nb_motifs),tens_Conf(:,:,:,nb_motifs));
