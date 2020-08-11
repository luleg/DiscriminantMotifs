function [ p,r ] = Amalgamation( A,p,r)
gamma = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function looks for improving the quality measure of a given
% clusterization in amlgaming some clusters.
%
% Inputs :
% A  : the matrix we are clusterizing
% p  : the reordered indices of the rows/columns of A, to make the clusters
%   appear
% r  : indicates the numbers of clusters [length(r)-1], and the size of
%   these clusters [for the cluster i : Cardinal(i) = r(i+1)-r(i)]
%
% Outputs :
% p and r for the new clusterisation
%
mn = size(A,1);

deux_m = sum(sum(A));

% Kc : number of clusters
Kc = length(r)-1;
p = p(:)'; r = r(:)';

% If we are working on the A rows
V = zeros(mn,Kc);
% V is a 'mn'x'Kc' matrix. Each column i of V indicates indices of the
% elements in the ith Cluster : V(k,i) = 1 if the kth element is in
% the ith cluster, and 0 if it isn't
for k = 1:Kc
    V(r(k):r(k+1)-1, k) = 1;
end
M = sum(A);
M = M(:);
% Sc : Condensed matrix : practical to work iteratively
A = A(p,p)';

Sc = 1/deux_m*V'*(A-gamma/deux_m*(M*M'))*V;

% Quality measure associated to the initial clusterisation
Qmes = trace(Sc);
clear V M

% The quality measure if two clusters k and k' are merged is given by :
% Q(k,k') = Qinit+2/mn(Sc(k,k')-Card(k,k')). So, we create the matrix
% G = 2/mn(Sc(k,k')-Card(k,k')) and while there is a non diagonal entry
% which is positive, we actually merge the clusters corresponding to its
% indices
G = Sc; G = G-diag(diag(G));
% while we can improve the quality measure, we keep on amalgaming
Continuer = true;
while Continuer
    %  we look for the maximal entry of G
    [Qref,i]=max(G);
    [Qref,j]=max(Qref);
    i=i(j);
    % we want to have i<j, in order to simplify the merge in terms of p and
    % r. Because the matrix G is symmetric, it changes nothing to exchange
    % i and j.
    if(i>j)
        aux=j; j=i;i=aux;
        % if the two indices are equal, it means that the max of G is 0. It
        % means that merging two clusters will not improve the quality measure,
        % so we stop the process. This test is done here to avoid computation
        % error
    elseif(i==j)
        Qref=-1; i=-1;j=-1;
    end
    % If none of the couple can improve the quality measure, we leave the
    % loop and so the function
    if(Qref <= 0)
        Continuer = false;
        % But if we have a couple which improves the quality measure, we
        % amalgame the two initial clusters and we update the variables
    else
        % We have find a couple (i,j) to amalgame, so we put the elements of
        % the jth cluster in the ith cluster, and we remove the jth cluster.
        Kc = Kc-1;
        % Here is the new quality measure (after amalgamation)
        Qmes = Qmes + 2*Qref;
        % We update the condensed matrix : we add jth column and row to the
        % ith ones and remove them.
        G(i,:)=G(i,:)+G(j,:);
        G(:,i)=G(i,:)';
        G(j,:)=[];G(:,j)=[];
        G = G -diag(diag(G));
        % we move the indices of the elements of the new cluster in a same block
        % this can be done since we always have i<j  (comblnk ensures this)
        p = [p(1:r(i+1)-1), p(r(j):r(j+1)-1), p(r(i+1):r(j)-1), p(r(j+1):r(end)-1)];
        % we finally update the r vector, to remove jth cluster, and change
        % the number of elements of ith cluster.
        nbelts = r(j+1)-r(j); r(j) = []; r(i+1:j-1) = r(i+1:j-1)+nbelts;
    end
end
p = p(:); r = r(:);

end
