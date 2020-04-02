% script to show the dispersion around classical motifs of networks from
% each field.
% Please note that this code is not as cleaned up as the others
clear variables;clc;
close all

% some global variables


%% Processing of networks
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

nb_reseaux = nb_fw+nb_elec+nb_stac+nb_soc;
fw  = 1:nb_fw;
elec= 1:nb_elec;elec= elec+nb_fw;
stac= 1:nb_stac;stac= stac+nb_elec+nb_fw;
soc = 1:nb_soc; soc = soc +nb_stac+nb_elec+nb_fw;

lreseaux = {};

for i=1:nb_fw
    lreseaux{i} = lreseaux_fw{i};
end
for i=1:nb_elec
    lreseaux{i+nb_fw} = lreseaux_elec{i};
end
for i=1:nb_stac
    lreseaux{i+nb_fw+nb_elec} = lreseaux_stac{i};
end
for i=1:nb_soc
    lreseaux{i+nb_fw+nb_elec+nb_stac} = lreseaux_soc{i};
end

V = generate_basis(lreseaux);


%% Processing of motifs. 
% You should uncomment the motif that interests you. You can add other kind
% of motifs you may be interested in, by adding its number of nodes (3/4),
% and its id as explained in Sec3 from Supplementary Material
% The formal process is
% name_motif = id_motifk == id_motif; % (k = 3 or 4)
% name_motif = indk(name_motif); % (k = 3 or 4)
% keep_motifs = [keep_motifs;name_motif+t]; % (t=0 if k==3; t=13 if k==4)

load('id_motifs.mat');
keep_motifs = [];
ind3 = 1:13;ind4 = 1:199;

% % Well-known motifs:

% % two bi-dir : 3-78
% bidir = id_motif3 == 78;
% bidir = ind3(bidir);
% keep_motifs = [keep_motifs;bidir];

% % bifan : 4-204
% bifan = id_motif4 == 204;
% bifan = ind4(bifan);
% keep_motifs = [keep_motifs;13+bifan];
 
% % trunc-bifan : 4-76
% bifan = id_motif4 == 76;
% bifan = ind4(bifan);
% keep_motifs = [keep_motifs;13+bifan];

% % feedforward loop : 3-38
ffl = id_motif3 == 38;
ffl = ind3(ffl);
keep_motifs = [keep_motifs;ffl];

% % long omnivory : 4-472
% lffl = id_motif3 == 12;
% lffl = ind3(lffl);
% keep_motifs = [keep_motifs;lffl];
 
% % bi-parallel : 4-904
% bipar = id_motif4 == 904;
% bipar = ind4(bipar);
% keep_motifs = [keep_motifs;13+bipar];
 
% % short feedback loop : 3-98
% sfbl = id_motif3 == 98;
% sfbl = ind3(sfbl);
% keep_motifs = [keep_motifs;sfbl];
 
% % complete triad : 3-238
% comp_tri = id_motif3 == 238;
% comp_tri = ind3(comp_tri);
% keep_motifs = [keep_motifs;comp_tri];
 
% % long feedback loop : 4-4740
% lfbl = id_motif4 == 4740;
% lfbl = ind4(lfbl);
% keep_motifs = [keep_motifs;13+lfbl];


id_motifs = [id_motif3;id_motif4];
indic_motifs = (keep_motifs>13)*4+(keep_motifs<14)*3;

indic_motifs = [indic_motifs,id_motifs(keep_motifs)];


nb_motifs = length(keep_motifs);
for i = 1:nb_motifs
    % figure(i) : distribution of networks in a boxplot view
    % figure(1000+i) : pointwise distribution of networks
    % figure(500+i) : graph corresponding to the motif id
    disp_dispersion(V(keep_motifs(i),:),fw,elec,stac,soc,indic_motifs(i,:),i)
end



function [] = disp_dispersion(V,fw,elec,stac,soc,motif,num_fig)
nb_fw = length(fw);nb_elec = length(elec);
nb_stac = length(stac);nb_soc = length(soc);
nb_reseaux = nb_fw+nb_elec+nb_stac+nb_soc;
vmeanfw = mean(V(fw));
vmeanelec = mean(V(elec));
vmeanstac = mean(V(stac));
vmeansoc = mean(V(soc));
vmean_f =mean([vmeanfw,vmeanelec,vmeanstac,vmeansoc]);
vmean_a = mean(V);


vmedfw = median(V(fw));
vmedelec = median(V(elec));
vmedstac = median(V(stac));
vmedsoc = median(V(soc));
vmedian_f =median([vmedfw,vmedelec,vmedstac,vmedsoc]);
vmedian_a = median(V);

indices = zeros(1,nb_reseaux);
indices(elec)=1;
indices(stac)=2;
indices(soc)=3;


figure(num_fig),clf
boxplot(V(fw),indices(fw),'whisker',0,'color','b','positions',0,'symbol','b+','widths',0.75)
hold on,
set(gca, 'box', 'on');
boxplot(V(elec),indices(elec),'whisker',0,'color','m','positions',1,'symbol','mo','widths',0.75)
boxplot(V(stac),indices(stac),'whisker',0,'color','k','positions',2,'symbol','kx','widths',0.75)
boxplot(V(soc),indices(soc),'whisker',0,'color','r','positions',3,'symbol','r*','widths',0.75),grid on

set(findobj(gca,'type','line'),'linew',1.5)
xlim([-1,4])
ylim([0 max(V)])

figure(num_fig), hold on,
plot([-1 4],[vmean_a vmean_a],'color',[0.8,0.5,0.5],'linestyle','--','linewidth',1.5)
figure(num_fig), hold on,
plot([-1 4],[vmean_f vmean_f],'color',[0.7,0.7,0.7],'linestyle','-.','linewidth',1.5)

if false
    figure(num_fig), hold on,
    plot([-1 4],[vmedian_a vmedian_a],'color',[0.8,0,0.5],'linestyle','--','linewidth',1.5)
    figure(num_fig), hold on,
    plot([-1 4],[vmedian_f vmedian_f],'color',[0.5,0,0.8],'linestyle','-.','linewidth',1.5)
end


if false
    ax= gca;
    ax.XTick = [0,1,2,3];
    ax.XTickLabel = {'Food Webs','Elec. Circ.','Disc. Struc','Social Net.'};
    ax.XTickLabelRotation = 75;
    ax.YLabel.String = 'Motif occurence number';
    ax.YLabel.FontSize = 17;
    lgd = legend('dataset mean','mean of field means','dataset median','median of field medians');lgd.FontSize = 15;
    ttl = title(['Disp. around ',int2str(motif(1)),'-',int2str(motif(2))]);
    ttl.FontSize = 17;
else
    ax = gca;
    ax.XTick = [];
end


%%

figure(1000+num_fig),clf
plot([vmean_a vmean_a],[-1 2],'color',[0.8,0.5,0.5],'linestyle','--','linewidth',1.5)
figure(1000+num_fig), hold on,
plot([vmean_f vmean_f],[-1,2],'color',[0.7,0.7,0.7],'linestyle','-.','linewidth',1.5)

if false
    figure(1000+num_fig), hold on,
    plot([vmedian_a vmedian_a],[-1 2],'color',[0.8,0,0.5],'linestyle','--','linewidth',1.5)
    figure(1000+num_fig), hold on,
    plot([vmedian_f vmedian_f],[-1 2],'color',[0.5,0,0.8],'linestyle','-.','linewidth',1.5)
end
hold on
plot(V(fw),2*ones(length(fw),1),'b*','markersize',7,'linewidth',1);
plot(V(elec),1*ones(length(elec),1),'m*','markersize',7,'linewidth',1);
plot(V(stac),0*ones(length(stac),1),'k*','markersize',7,'linewidth',1);
plot(V(soc),-1*ones(length(soc),1),'r*','markersize',7,'linewidth',1);


ax= gca;
ax.YTick = [-1,0,1,2];
ax.YTickLabel = {'Soc. Net.','Disc.Struc.','Elec. Circ.','Food Webs'};

% ax.XLabel.String = 'Motif occurence number';
% ax.XLabel.FontSize = 11;
% 
% ttl = title([int2str(motif(1)),'-',int2str(motif(2))]);
% ttl.FontSize = 20;
%%
print_graph(motif(1),motif(2),500+num_fig);

end

function [G] = print_graph(node,motif,num_fig)
id = motif;
if node == 3
    power_2 = [0 128 64 32 0 8 4 2 1];
elseif node == 4
    power_2 = [0 16384 8192 4096 2048 1024 512 256 128 64 0 16 8 4 2 1];
end

A = zeros(1,node*node);
cpt = 1;
while motif>0
    tmp = mod(motif,power_2(1));
    power_2(1) = [];
    if tmp < motif
        A(cpt) = 1;
        motif = tmp;
    end
    cpt = cpt+1;
end

A = fliplr(A);
A = reshape(A,[node node]);
A = A';
G = digraph(A);
if nargin <3, num_fig = 1000;end
figure(num_fig),clf
pl = plot(G);


if node == 3
    pl.XData = [0 -1 1];
    pl.YData = [1 0 0];
elseif node == 4
    pl.XData = [0 1 1 0];
    pl.YData = [1 1 0 0];
end
ttl = title(['id : ',int2str(id)]);
ttl.FontSize = 20;
pl.MarkerSize = 20;
pl.NodeFontSize = 30;
pl.NodeFontWeight = 'bold';

pl.LineWidth = 5;
pl.ArrowSize = 30;

end


