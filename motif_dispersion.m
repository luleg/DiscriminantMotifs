% script to show the dispersion around discriminant motifs of networks from
% each field.
clear variables;clc;
close all

% some global variables

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
V = V(keep_motifs,:);


load('id_motifs.mat');
id_motifs = [id_motif3;id_motif4];


indic_motifs = (keep_motifs>13)*4+(keep_motifs<14)*3;
indic_motifs = [indic_motifs, id_motifs(keep_motifs)];
%
% id_motifs = [id_motif3;id_motif4];
% indic_motifs = [[3*ones(13,1);4*ones(199,1)],id_motifs];
% nb_motifs = 212;
for i = 1:nb_motifs
    disp_dispersion(V(i,:),fw,elec,stac,soc,indic_motifs(i,:),i)
end


% disp_dispersion(V(9,:),fw,elec,stac,soc,indic_motifs(9,:),i)


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

figure(num_fig), hold on,
plot([-1 4],[vmedian_a vmedian_a],'color',[0.8,0,0.5],'linestyle','--','linewidth',1.5)
figure(num_fig), hold on,
plot([-1 4],[vmedian_f vmedian_f],'color',[0.5,0,0.8],'linestyle','-.','linewidth',1.5)



if (num_fig==1)
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
    ttl = title([int2str(motif(1)),'-',int2str(motif(2))]);
    ttl.FontSize = 20;
end


%%

figure(1000+num_fig),clf
plot([vmean_a vmean_a],[-1 2],'color',[0.8,0.5,0.5],'linestyle','--','linewidth',1.5)
figure(1000+num_fig), hold on,
plot([vmean_f vmean_f],[-1,2],'color',[0.7,0.7,0.7],'linestyle','-.','linewidth',1.5)

figure(1000+num_fig), hold on,
plot([vmedian_a vmedian_a],[-1 2],'color',[0.8,0,0.5],'linestyle','--','linewidth',1.5)
figure(1000+num_fig), hold on,
plot([vmedian_f vmedian_f],[-1 2],'color',[0.5,0,0.8],'linestyle','-.','linewidth',1.5)

hold on
plot(V(fw),2*ones(length(fw),1),'b*','markersize',7,'linewidth',1);
plot(V(elec),1*ones(length(elec),1),'m*','markersize',7,'linewidth',1);
plot(V(stac),0*ones(length(stac),1),'k*','markersize',7,'linewidth',1);
plot(V(soc),-1*ones(length(soc),1),'r*','markersize',7,'linewidth',1);
ax = gca;
ax.XTick = [];
ttl = title([int2str(motif(1)),'-',int2str(motif(2))]);
ttl.FontSize = 20;
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


