% Script that compute the PCA of basis V
% V = [v1,...,vn], where vk is the kth network characterised by its
% decomposition onto k-node motifs
[t_vect,nb_tr] = size(V);
% compute the mean element of the dataset V
x_b = 1/nb_reseaux*V*ones(nb_tr,1);
X = V-x_b*ones(1,nb_tr); % X is the centered dataset


% Variance-Covariance matrix of V:
S = 1/nb_reseaux * X*X';
% just a trick to enforce S to be symmetric, hence eig provides real
% values (S is expected to be symmetric to the rounding errors)
S = 1/2*(S+S');
% [U,D,flag] = eigs(S,20,'lr'); % U : basis of eigenvectors, D : diagonal matrix containing
[U,D,flag] = eigs(S,min(t_vect,20),'lr'); % U : basis of eigenvectors, D : diagonal matrix containing
if flag ~=0
    disp ("Imprecision in Eigenvectors")
end
% eigenvalues
[D,ind_sort_D] = sort(diag(D),'descend');
U = U(:,ind_sort_D);

% takes the indexes required to obtain the target percentage of trace (p)
if p>0
    tot = sum(D); cD = cumsum(D);
    ind_sort_D = 1:length(cD);
    ind_sort_D = ind_sort_D(cD>p*tot);
    nb_ax = max(ind_sort_D(1),3);
    prct_trace = nb_ax;
else
    prct_trace = sum(D(1:nb_ax))/sum(D);
end


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
    pause
    
end

% We compute the coordinates of the individuals in the new basis, with the
% limited number of components
C = X'*U(:,1:nb_ax);

if disp_fig
    disp_dataset(C,fw,elec,stac,soc,2);
end