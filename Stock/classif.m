k1=1;kf = 3;
[p,n] = size(C);
D=sqdist(C,C);
dmax = max(max(D));
sigma = dmax/n^(1/p);
sigma=sigma/8;
S = exp(-D/sigma^2);
S=1/2*(S+S');
S=S-diag(diag(S));
D = S*ones(n,1);
L=diag(1./sqrt(D))*S*diag(1./sqrt(D));
L=1/2*(L+L');
[U,D] = eig(L);
[d,inds]=sort(diag(D),'descend');
U=U(:,inds);
X=U(:,k1:kf);
nX = sqrt((X.*X)*ones(kf-k1+1,1));

Y=diag(1./nX)*X;