% Script that compute the PCA of basis V
% V = [v1,...,vn], where vk is the kth network characterised by its
% decomposition onto 3-node and 4-node motifs

% compute the mean element of the dataset V
x_b = 1/nb_reseaux*V*ones(nb_reseaux,1);
X = V-x_b*ones(1,nb_reseaux); % X is the centered dataset

% Display the motif dispersion of networks from each field around the mean
% element
if disp_fig
    load('id_motifs.mat')
    mot3_om = 1:13;
    labx = id_motif3;
    mot4_om = 1:199;limbar3 = [0.5,13.5];
    labx4 = id_motif4;limbar4 = [0.5,199.5];
    X3 = X(1:13,:);X4 = X(14:end,:);
    
    figure(100),clf
    sgtitle('Dispersion around mean element for Foodwebs')
    subplot(2,1,1),plot([1 13],[0 0],'k-.'),hold on,plot(X3(mot3_om,fw)),xlim([1,13]),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),plot([1 199],[0 0],'k-.'),hold on,plot(X(14:end,fw)),xlim([1,199]),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
    
    figure(110),clf
    sgtitle('Dispersion around mean element for the mean of foodweb')
    subplot(2,1,1),bar(mean(X3(mot3_om,fw),2)),xlim(limbar3),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),bar(mean(X(14:end,fw),2)),xlim(limbar4),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')

    
    figure(101),clf
    sgtitle('Dispersion around mean element for Electronic Circuits')
    subplot(2,1,1),plot([1 13],[0 0],'k-.'),hold on,plot(X3(mot3_om,elec)),xlim([1 13]),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),plot([1 199],[0 0],'k-.'),hold on,plot(X(14:end,elec)),xlim([1 199]),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
    
    figure(111),clf
    sgtitle('Dispersion around mean element for the mean of electronic circuit')
    subplot(2,1,1),bar(mean(X3(mot3_om,elec),2)),xlim(limbar3),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),bar(mean(X(14:end,elec),2)),xlim(limbar4),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
    
    
    figure(102),clf
    sgtitle('Dispersion around mean element for Discourse Structures')
    subplot(2,1,1),plot([1 13],[0 0],'k-.'),hold on,plot(X3(mot3_om,stac)),xlim([1 13]),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),plot([1 199],[0 0],'k-.'),hold on,plot(X(14:end,stac)),xlim([1 199]),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
    
    figure(112),clf
    sgtitle('Dispersion around mean element for the mean of discourse structure')
    subplot(2,1,1),bar(mean(X3(mot3_om,stac),2)),xlim(limbar3),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),bar(mean(X(14:end,stac),2)),xlim(limbar4),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
    
    
    figure(103),clf
    sgtitle('Dispersion around mean element for Social Networks')
    subplot(2,1,1),plot([1 13],[0 0],'k-.'),hold on,plot(X3(mot3_om,soc)),xlim([1 13]),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),plot([1 199],[0 0],'k-.'),hold on,plot(X(14:end,soc)),xlim([1 199]),xticks(mot4_om);xticklabels(labx4);
    
    figure(113),clf
    sgtitle('Dispersion around mean element for the mean of social network')
    subplot(2,1,1),bar(mean(X3(mot3_om,soc),2)),xlim(limbar3),xticks(mot3_om);xticklabels(labx);
    title('For 3-node motifs')
    subplot(2,1,2),bar(mean(X(14:end,soc),2)),xlim(limbar4),xticks(mot4_om);xticklabels(labx4);
    title('For 4-node motifs')
end

% Variance-Covariance matrix of V:
S = 1/nb_reseaux * X*X';

% just a trick to enforce S to be symmetric, hence eig provides real
% values (S is expected to be symmetric to the rounduing errors)
S = 1/2*(S+S');
[U,D] = eig(S); % U : basis of eigenvectors, D : diagonal matrix containing
% eigenvalues

[D,ind_sort_D] = sort(diag(D),'descend');
U = U(:,ind_sort_D);

% takes the indexes required to obtain the target percentage of trace (p)
tot = sum(D); cD = cumsum(D);
ind_sort_D = 1:length(cD);
ind_sort_D = ind_sort_D(cD>p*tot);
nb_ax = max(ind_sort_D(1),3); % We keep at least 3 principal axes (to
% be able to plot the 3D figures)

% Display the eigenvalues in descending order, and the cumulative sum, and
% the threshold that shows the number of eigenvalues kept
if disp_fig
    figure(1),clf,
    sgtitle('Percentage of trace')
    subplot(2,1,1)
    plot(D/sum(D)*100,'r-*')
    xlim([1,min(t_vect,nb_reseaux)])
    title('in descending order')
    subplot(2,1,2)
    plot(cumsum(D)/sum(D)*100,'r-*')
    xlim([1,min(t_vect,nb_reseaux)])
    title('Cumulative sum')
    figure(1),hold on,
    subplot(2,1,1),hold on,
    plot([nb_ax+0.5 nb_ax+0.5],[0 max(D/sum(D)*100)],'k-.','linewidth',2)
    legend('eigenvalues','eigenvalues kept')
    subplot(2,1,2),hold on,
    plot([nb_ax+0.5 nb_ax+0.5],[min(cumsum(D)/sum(D)*100) 100],'k-.','linewidth',2)
    
end

% We compute the coordinates of the individuals in the new basis, with the
% limited number of components
C = X'*U(:,1:nb_ax);

if disp_fig
    disp_dataset(C,fw,elec,stac,soc,2);
end