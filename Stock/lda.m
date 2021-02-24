% Script that compute the LDA of basis V
% V = [v1,...,vn], where vk is the kth network characterised by its
% decomposition onto k-node motifs
nb_tr = size(V,2);
% compute the mean element of the dataset V...
x_b = 1/nb_tr*V*ones(nb_tr,1);
% ... and of each benchmark.
mu_fw = 1/nb_basis*V(:,fw)*ones(nb_basis,1);
mu_elec = 1/nb_basis*V(:,elec)*ones(nb_basis,1);
mu_soc = 1/nb_basis*V(:,soc)*ones(nb_basis,1);
mu_stac = 1/nb_basis*V(:,stac)*ones(nb_basis,1);
% Between Scatter Matrix:
Btmp = sqrt(nb_basis)*([mu_fw,mu_elec,mu_stac,mu_soc]-x_b);
Sb = Btmp*Btmp';
Sb = 1/2*(Sb+Sb');

% Within Scatter Matrix:
mmean = ones(size(V));
ee = ones(1,nb_basis);
mmean(:,fw) = mu_fw*ee;
mmean(:,elec) = mu_elec*ee;
mmean(:,stac) = mu_stac*ee;
mmean(:,soc) = mu_soc*ee; 
Wtmp = V-mmean;
Sw=Wtmp*Wtmp';
Sw=1/2*(Sw+Sw');

% Solutions of the generalized eigenvalue problem :
[U,D] = eig(pinv(Sw)*Sb);
[d,inds] = sort(real(diag(D)),'descend');
U = U(:,inds);
U = U(:,1:nb_ax);

% Gram-schmidt to have an orthonormal basis for projection.
Q = zeros(size(U));
for k=1:size(U,2)
    Q(:,k)=U(:,k);
    for j=1:k-1
        Q(:,k)=Q(:,k)-(U(:,k)'*Q(:,j))*Q(:,j);
    end
    Q(:,k) = Q(:,k)/norm(Q(:,k));
end
U = Q;

C = (V-x_b)'*U;
