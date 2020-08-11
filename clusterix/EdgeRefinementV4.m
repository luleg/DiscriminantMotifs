function [ p,r,m1 ] = EdgeRefinementV4( U,r,nb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDGEREFINEMENTV4 is a process that corrects a 
% clustering proposed by a convolution analysis by
% moving some edges returned by the filter to their 
% exact correct position in the vector (may be in a 
% neighbourhood), potentially sharpens the vector by
% means of projection onto the basis of several 
% eigenvector that may convey additional
% information about the separation between two 
% clusters, and check if the found edges are sharp 
% enough.
% INPUT
%    U	: a basis of some leading eigenvectors
%    r  : the separation found by the filter
%    nb	: the indice of the analysed vector as column
%  of U
% OUTOUT
%    p,r : the clustering (p(r(i):r(i+1)-1) are the
% indices of the i-th block
%    m1  : the confidence we have in the returned 
% clustering (negative if we are strongly confident,
%  equal to the measure of sharpness of the best edge 
% otherwise)


%% 3 for threshold, 4 for knn
tol = 3; % to filter spurious edges (> 2) (the level of confidence we want to have)
u = U(:,nb); [~,p] = sort(u); 
kc = length(r)-1;
% if no separation found by the filter
if length(r)==2
    m1 =0;
    return
end

% create the basis of characteristic vector : V(i,j) != 0 means that 
% the i-th element of u belongs to the j-th block. 
V = zeros(r(end)-1,kc);
for k=1:kc
    V(r(k):r(k+1)-1,k) = 1/sqrt(r(k+1)-r(k));
end
V(p,:) = V;

% vk is a staircase vector expected to be close to u
vk = V*(V'*u);


for j=1:5
    for i=1:4
        % We update the clustering by finding the staircase vector that
	% fits the best with u
        [r,p] = MAJ(u,vk);
        
        % we then create the corresponfing basis of characteristic vectors
        kc = length(r)-1;
        V = zeros(r(end)-1,kc);
        for k=1:kc
            V(r(k):r(k+1)-1,k) = 1/sqrt(r(k+1)-r(k));
        end
        V(p,:) = V;
        vk = V*(V'*u);

    end
    if(kc>1)
        % we look for the vector in the basis U that fits the best with the 
	% staircase vector that represents the clustering. 
        u = U*(U'*vk); u = u /norm(u);
        vk = V*(V'*u);
    else
        break
    end
    
end


% filter the spurious edges :
ratio = sharpness(u,r); % computation of their sharpness measure
ratio = ratio(:)';

non_spur = ones(1,length(ratio)+2);
non_spur(2:end-1) = (ratio>=tol);

% we look for the confidence we have in the edge(s) we return
[m1,ind_max] = max(ratio);
if(m1>=tol),m1=-2; % we are confident
elseif(isempty(m1)),m1=0;% we are not confident
end
non_spur(ind_max+1)=true;

r = r(logical(non_spur));

end

%% Sharpness measure
function [ratio] = sharpness(u,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHARPNESS measure the sharpness measure of a set of edges
% provided by r on a signal u.

kc = length(r)-1;
u = sort(u);
alpha = zeros(1,kc);err = zeros(1,kc);
for k=1:kc
    alpha(k) = mean(u(r(k):r(k+1)-1));
    err(k)   = 1/(r(k+1)-r(k))*sum(abs(u(r(k):r(k+1)-1)-alpha(k)));
end
ratio = (alpha(2:end)-alpha(1:end-1))./(err(2:end)+err(1:end-1));
end


%% Clustering update
function [r,p] = MAJ(u,v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAJ updates the clustering provided by a staircase vector 
% v according to the values of u.
[u,p] = sort(u(:));

v = unique(v);
val_seuil = (v(2:end)+v(1:end-1))/2;
val_seuil = [val_seuil(:)',val_seuil(end)+10];
r = [1]; k = 1;
seuil = val_seuil(k);

for i = 1:length(u)
    if(u(i)>seuil)
        r = [r i];
        while(u(i)>seuil)
            k = k+1;
            seuil = val_seuil(k);
        end
    end
end
r = [r i+1];
end
