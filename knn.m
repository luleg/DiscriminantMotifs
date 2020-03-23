% knn 

% creation of the labels of the learning basis
class_learn = ones(1,nb_reseaux);
class_learn(fw) = 1;
class_learn(elec) = 2;
class_learn(stac) = 3;
class_learn(soc) = 4;

% attribution of the field to networks from the test basis
classek = fct_knn(C',Ct',class_learn);

% display of the attribution of the test networks from each field
if disp_fig
    figure(100+nb_k),clf,plot(classek,'*','linewidth',2)
    title(['Labels attributed by the ',int2str(nb_k),'-nn classifier'])
    hold on, plot([nbt_fw+1/2 nbt_fw+1/2],[1 4],'k-.',...
        [nbt_fw+nbt_elec+1/2 nbt_fw+nbt_elec+1/2],[1 4],'k-.',...
        [nbt_fw+nbt_elec+nbt_stac+1/2 nbt_fw+nbt_elec+nbt_stac+1/2],[1 4],'k-.','linewidth',2)
    ax = gca;    
    ax.YTick = [1,2,3,4];
    ax.YTickLabel = {'Food Webs','Elec. Circ','Disc. Struct','Social Net.'};
    ax.YLim = [0.9,4.1];
    ax.YLabel.String= 'attributed label';
    ax.XTick = [mean(fwt),mean(elect),mean(stact),mean(soct)];
    ax.XTickLabel = {'Food Webs','Elec. Circ','Disc. Struct','Social Net.'};
    ax.XLim = [1, nbt_reseaux];
    ax.XLabel.String = 'true label';
end

[Conf,infos] = compute_scores(classek,fwt,elect,stact,soct);