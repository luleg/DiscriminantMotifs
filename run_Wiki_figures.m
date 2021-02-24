close all;
clear;
clc

global p
p =-1; 
global nb_ax
nb_ax =4;
global disp_fig
disp_fig = 0; % do you want to display the Figures ?

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
V_quali = generate_basis(lreseaux,1,0,0);

lreseaux = cell(1,nb_nn);
for i=1:nb_nn
    lreseaux{i} = lreseaux_nn{i};
end
V_nn= generate_basis(lreseaux,1,0,0);

lreseaux = cell(1,nb_ge);
for i=1:nb_ge
    lreseaux{i} = lreseaux_ge{i};
end
V_ge= generate_basis(lreseaux,1,0,0);


V_pp = [V_nn,V_ge];
nbBasis = 20;

V = [V_quali,V_pp];
pca


axes = [1,3];

ind_quali = 1:nb_quali;
ind_nn = 1:nb_nn;ind_nn= ind_nn+nb_quali;
ind_ge = 1:nb_ge; ind_ge = ind_ge+nb_quali+nb_nn;
ind_pp = [ind_nn,ind_ge];

figure(1),clf, hold on,grid on
plot(C(ind_quali,axes(1)),C(ind_quali,axes(2)),'ms','markersize',10,'linewidth',1.5)
plot(C(ind_pp,axes(1)),C(ind_pp,axes(2)),'bo','markersize',10,'linewidth',1.5)


P = polyfit(C(ind_quali,1),C(ind_quali,3),2);

xx = -0.4:0.01:0.3;
yy = P(1)*xx.^2+P(2)*xx+P(3);

figure(1),hold on, plot(xx,yy,'m:','linewidth',2)

P = polyfit(C(ind_pp,1),C(ind_pp,3),2);

xx = -0.4:0.01:0.6;
yy = P(1)*xx.^2+P(2)*xx+P(3);

figure(1),hold on, plot(xx,yy,'b:','linewidth',2)

inds = 1:88;
ii = inds(V(1,:)>0.985);
figure(1),hold on,
plot(C(ii,axes(1)),C(ii,axes(2)),'k+','linewidth',1.5)
ax = gca;
ax.XLabel.String = '1st ppal axis';
ax.XLabel.FontSize=15;
ax.XLabel.FontWeight='bold';
ax.YLabel.String = '3rd ppal axis';
ax.YLabel.FontSize=15;
ax.YLabel.FontWeight='bold';
lgd = legend('quality','contested','fit quality','fit contested');
lgd.FontSize = 10;
lgd.FontWeight = 'bold';


%%
nbMot=13;
load('Stock/id_motifs.mat');
indic_motif = [3*ones(nbMot,1),id_motif3];


figure(2),clf
plt = plot(V_quali);
ax = gca;
ax.XAxis.TickValues = 1:nbMot;
ind_mot_str = {};
for k=1:nbMot
    ind_mot_str{k} = [int2str(indic_motif(k,1)),'-',int2str(indic_motif(k,2))];
end
ax.XAxis.TickLabels  = ind_mot_str;
ax.XAxis.TickLabelRotation = 70;
ax.XAxis.FontSize = 15;
xlim([1 13])
hold on, grid on

ax.YAxis.Label.String = 'Embeddings';
ax.YAxis.Label.FontSize= 20;
ax.Title.String = 'Quality Articles';
ax.Title.FontSize=25;
ax.Title.FontWeight='bold';


figure(3),clf
plt = plot(V_pp);
ax = gca;
ax.XAxis.TickValues = 1:nbMot;
ind_mot_str = {};
for k=1:nbMot
    ind_mot_str{k} = [int2str(indic_motif(k,1)),'-',int2str(indic_motif(k,2))];
end
ax.XAxis.TickLabels  = ind_mot_str;
ax.XAxis.TickLabelRotation = 70;
ax.XAxis.FontSize = 15;
xlim([1 13])
hold on, grid on

ax.YAxis.Label.String = 'Embeddings';
ax.YAxis.Label.FontSize= 20;
ax.Title.String = 'Contested Articles';
ax.Title.FontSize=25;
ax.Title.FontWeight='bold';

%%
edges = 0:0.025:1;
quali_14 = V_quali(3,:);
pp_14 = V_pp(3,:);

figure(4),clf
histogram(quali_14,edges,'normalization','cdf')
hold on,histogram(pp_14,edges,'normalization','cdf')
ttl = title("CDF of Graphlet 3-14");
ttl.FontSize = 25;
lgd = legend("Quality Articles","Contested Articles");
lgd.FontSize =15;
lgd.FontWeight='bold';
ax = gca;
ax.XLim = [0 0.7];
ax.YLim = [0 1];
ax.XLabel.String = 'Proportion of 3-14 graphlet';
ax.XLabel.FontSize= 15;
ax.XLabel.FontWeight = 'bold';

ax.YLabel.String = 'Cummulative Density';
ax.YLabel.FontSize= 15;
ax.YLabel.FontWeight = 'bold';

edges = 0:0.025:1;
quali_74 = V_quali(7,:);
pp_74 = V_pp(7,:);

figure(5),clf
histogram(quali_74,edges,'normalization','cdf')
hold on,histogram(pp_74,edges,'normalization','cdf')
ttl = title("CDF of Graphlet 3-74");
ttl.FontSize = 25;
lgd = legend("Quality Articles","Contested Articles");
lgd.FontSize =15;
lgd.FontWeight='bold';
ax = gca;
ax.XLim = [0 0.5];
ax.YLim = [0 1];
ax.XLabel.String = 'Proportion of 3-74 graphlet';
ax.XLabel.FontSize= 15;
ax.XLabel.FontWeight = 'bold';

ax.YLabel.String = 'Cummulative Density';
ax.YLabel.FontSize= 15;
ax.YLabel.FontWeight = 'bold';