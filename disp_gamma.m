function [] = disp_gamma(Gamt,num_fig,threshold,disp_all_motifs,nb_mot,max_mot)
% Plot the Gamma Score in descending order
% Compulsory parameters:
% Gamt : the gamma score of all motifs
% num_fig: identifier of the figure
% Optional parameters:
% threshold : number of principal axes used to compute the gamma score (to 
%           plot in the figure title). Default : title without this info
% disp_all_motifs : if true, all motif identifiers are given on the x axis
%           Default : false
% nb_mot : number of motifs with the ihghest gamma score whose label should
%           be displayed with tick. Default : 10
% max_mot : right limit of the x axis (only the max_mot motifs with highest
%           gamma score are displayed). Default : 20


plot_title = true;

switch(nargin)
    case 2
        plot_title = false;
        disp_all_motifs = false;
        nb_mot = 10;
        max_mot = 20;
    case 3
        disp_all_motifs = false;
        nb_mot = 10;
        max_mot = 20;
    case 4
        nb_mot = 10;
        max_mot = 20;
    case 5
        max_mot = 20;
end

if threshold == -1 
    plot_title = false;
end
load('id_motifs.mat');
id_motifs = [id_motif3;id_motif4];

[Gamt,perm_mot] = sort(Gamt,'descend');

perm_mot = perm_mot(1:max_mot);
indic_motif = (perm_mot>13)*4+(perm_mot<14)*3;
indic_motif = [indic_motif, id_motifs(perm_mot)];

figure(num_fig),clf,
plt = plot(Gamt,'r-*','linewidth',2);grid on
xlim([1 max_mot])
for i =1:nb_mot
    figure(num_fig), hold on,
    plt = plot(i,Gamt(i),'r*');
    datatip(plt,'Fontsize',14,'DataIndex',i);
    mlabel = [int2str(indic_motif((i),1)),'-motif ',int2str(indic_motif((i),2))];
    plt.DataTipTemplate.DataTipRows(1) = mlabel;
    plt.DataTipTemplate.DataTipRows(2) = '';
end

if plot_title
    titre = ['\gamma-score of motifs within the ',int2str(threshold),' principal axes'];
else
    titre = '\gamma-score of motifs';
    
end

figure(num_fig),hold on,
ttl = title(titre);
ttl.FontSize  = 15;

ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Label.String = '\gamma(motif_k)';
ax.YAxis.Label.FontSize = 15;

if disp_all_motifs
    ax.XAxis.TickValues = 1:max_mot;
    ind_mot_str = {};
    for k=1:max_mot
        ind_mot_str{k} = [int2str(indic_motif(k,1)),'-',int2str(indic_motif(k,2))];
    end
    ax.XAxis.TickLabels  = ind_mot_str;
    ax.XAxis.TickLabelRotation = 70;
    ax.XAxis.FontSize = 7;
else
    ax.XAxis.Label.String = 'k';
    ax.XAxis.Label.FontSize = 15;
end

end

