% script to discover the discriminant motifs based on the PCA of a subset
% of the networks
clear variables;clc;
close all

% some global variables

global nb_ites
nb_ites = 500; % number of iterations (number of classifiers based on 1 PCA
% basis instance)
global disp_fig
disp_fig = false; % do you want to display the Figures ?

do_save = true; % If true, saves the order of motifs according to the mean 
% gamma score, in a mat file called discriminant_motifs.mat

global p
p = 0.75; % percentage of trace to keep in the PCA

nb_basis = 40; % number of networks from each field used to generate the
% PCA basis
% (/!\ the smallest database is the electronic circuits : 52 networks)

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

global t_vect
t_vect = 199+13; % size of the vectors used to represent networks :
% 199 4-node motifs, 13 3-node motifs

tens_gamma = zeros(t_vect,nb_ites); % tens_gamma(i,k) will contain the
% gamma score of ith motif in the PCA provided by the kth iteration

tens_k = zeros(1,nb_ites); % number of pca to obatin the target % of trace


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

for iter = 1:nb_ites
    if mod(iter,50) == 0
        iter
    end
    
    % random generation of the basis
    ind_fw = randperm(nb_fwi);
    ind_elec = randperm(nb_eleci);
    ind_stac = randperm(nb_staci);
    ind_soc = randperm(nb_soci);
    
    %% PCA
    % Generation of the basis
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
    
    % calling of the acp script, that manages the PCA process
    pca
    
    tens_k(iter)=nb_ax; % Knowing how many component are klept at each
    %% Discover motifs
    discriminant_motifs
    tens_gamma(:,iter) = Gamt;
    
end

mgamma = mean(tens_gamma,2);

disp_gamma(mgamma,1,-1,true,10,30);

[mgamma,perm_mot] =sort(mgamma,'descend');
sort_motif = [perm_mot,mgamma]; 
if do_save
    save('discriminant_motifs','sort_motif')
end