function [] = disp_gamma_errorbar(Gamt,num_fig,max_mot,disp_all_motifs,nb_mot,threshold,ttl)
% Plot the Gamma Score in descending order
% Compulsory parameters:
% Gamt : the gamma score of all motifs : 
% Gamt is a matrix nb_ites x nb_motifs (212)
% 
% Optional parameters:
% num_fig: identifier of the figure
%           default : 1
% max_mot : right limit of the x axis (only the max_mot motifs with highest
%           gamma score are displayed). Default : 20
% disp_all_motifs : if true, all motif identifiers are given on the x axis
%           Default : true
% nb_mot : number of motifs with the ihghest gamma score whose label should
%           be displayed with tick. Default : 0
% threshold : number of principal axes used to compute the gamma score (to 
%           plot in the figure title). Default : title without this info


plot_title = true;

switch(nargin)
    case 1
        num_fig = 51;
        max_mot = 20;
        disp_all_motifs = true;
        nb_mot = 0;
        threshold = -1;
    case 2
        max_mot = 20;
        disp_all_motifs = true;
        nb_mot = 0;
        threshold = -1;
    case 3
        disp_all_motifs = true;
        nb_mot = 0;
        threshold = -1;
    case 4
        nb_mot = 0;
        threshold = -1;
    case 5
        threshold = -1;
    case 7
        threshold = -2;
end


load('id_motifs.mat');
id_motifs = [id_motif3;id_motif4];

mean_gamma = mean(Gamt);
[mean_gamma,perm_mot] = sort(mean_gamma,'descend');
perm_mot = perm_mot(:);
Gamt = Gamt(:,perm_mot);
std_gamma = std(Gamt);
Gamt = mean_gamma;

perm_mot = perm_mot(1:max_mot);
indic_motif = (perm_mot>13)*4+(perm_mot<14)*3;

indic_motif = [indic_motif, id_motifs(perm_mot)];

figure(num_fig),clf,
plt = errorbar(1:max_mot,Gamt(1:max_mot),std_gamma(1:max_mot),'r-*','linewidth',2);grid on
xlim([1 max_mot])
for i =1:nb_mot
    figure(num_fig), hold on,
    plt = plot(i,Gamt(i),'r*');
    datatip(plt,'Fontsize',14,'DataIndex',i);
    mlabel = [int2str(indic_motif((i),1)),'-motif ',int2str(indic_motif((i),2))];
    plt.DataTipTemplate.DataTipRows(1) = mlabel;
    plt.DataTipTemplate.DataTipRows(2) = '';
end

switch threshold
    case -1
        titre = '\gamma-score of graphlets';
    case -2
        titre = ttl;
    otherwise
        titre = ['\gamma-score of graphlets for ',int2str(threshold),' principal axes'];
end
        

figure(num_fig),hold on,
ttl = title(titre);
ttl.FontSize  = 20;

ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Label.String = '\gamma(graphlet_k)';
ax.YAxis.Label.FontSize = 20;

if disp_all_motifs
    ax.XAxis.TickValues = 1:max_mot;
    ind_mot_str = {};
    for k=1:max_mot
        ind_mot_str{k} = [int2str(indic_motif(k,1)),'-',int2str(indic_motif(k,2))];
    end
    ax.XAxis.TickLabels  = ind_mot_str;
    ax.XAxis.TickLabelRotation = 70;
    ax.XAxis.FontSize = 15;
else
    ax.XAxis.Label.String = 'k';
    ax.XAxis.Label.FontSize = 15;
end

end

