% script to discover the discriminant motifs based on the PCA of a subset
% of the networks

addpath('Stock');
clear variables;clc;
close all

load RandPermDatasets
nb_ites = size(lind_elec,1);
% nb_ites = 10;

% some global variables

global disp_fig
disp_fig = 0; % do you want to display the Figures ?

global do_save
do_save = true; % do you want to save the gamma-values of graphlets in a 
% mat file?

global p
p =-1; % percentage of trace to keep in the PCA (if <0 one keeps a number 
% of axes instead of a percentage of trace)
global nb_ax 
nb_ax = 11;

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
norm_2 = 1; % second step (normalisation in norm-2)
% 1 : x/norm(x)
% 2 : x(motifk)/norm([xi(motifk), i = 1 :nb_graphs])
% else :  no normalisation



tens_per = zeros(1,nb_ites); % number of ppal axes to obtain the target % 
% of trace (or % of trce given a fie number of ppal axes)
tens_Conf = zeros(4,4,nb_ites);
% tens_Conf(i,j,k) : number of networks from field i that has been
% classified in field j by the classifier at the kth iteration

global t_vect
t_vect = 199+13; % size of the vectors used to represent networks :
% 199 4-node motifs, 13 3-node motifs
tens_gamma = zeros(nb_ites,t_vect);

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


for iter = 1:nb_ites
    if mod(iter,50) == 0
        disp(['iteration number: ',int2str(iter),' over ',int2str(nb_ites),'.']);
    end
    
    % random generation of the basis
    ind_fw = lind_fw(iter,1:nb_basis);
    ind_elec = lind_elec(iter,1:nb_basis);
    ind_stac = lind_stac(iter,1:nb_basis);
    ind_soc = lind_soc(iter,1:nb_basis);
    
    %% PCA
    V = [V_fw(:,ind_fw(1:nb_basis)),V_elec(:,ind_elec(1:nb_basis)),V_stac(:,ind_stac(1:nb_basis)),V_soc(:,ind_soc(1:nb_basis))];
    
    % calling of the acp script, that manages the PCA process
    pca
    tens_per(iter)= prct_trace;
    
    %% Discover motifs
    discriminant_motifs_gamma2
    tens_gamma(iter,:) = Gamt;
    
end

if do_save
    save('Stock/gamma_to_analyse','tens_gamma','tens_per')
end

disp_gamma_errorbar(tens_gamma,51,20,true,0,nb_ax)