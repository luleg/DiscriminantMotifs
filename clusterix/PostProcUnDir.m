function [p,r] = PostProcUnDir(A,p,r,nb_bloc_preproc)
% POSTPROCUNDIR is the postprocessing of CLUSTERIX for undirected graph.
% It merges small clusters detected in the preprocessing stage with those
% detecting with spectral tools, by improving the modularity measure.
% INPUT
%    A	: the doubly-stochastic scaling of the adjacency matrix of the 
% (sub)graph we want to clusterize.
%    p,r : the clustering (p(r(i):r(i+1)-1) are the indices of the elements
% that belong to the i-th cluster)
%    nb_bloc_preproc : the number of cluster taken apart during the 
% preprocessing stage
% OUTPUT 
%    p,r : the clustering

[~,n] = size(A);
m = n/2;

% Kc : number of clusters
Kc = length(r)-1;
la_diag = 0:Kc:Kc*(Kc-1);
la_diag = la_diag +[1:Kc];
p = p(:)'; r = r(:)';

V = zeros(n,Kc);
for i=1:Kc
    V(p(r(i):r(i+1)-1),i)=ones(r(i+1)-r(i),1);
end

R = sum(A,2);R = R*R';

Sc = V'*(A-R/(2*m))*V;Sc = Sc/(2*m);
% Quality measure associated to the initial clusterisation

Qmes = trace(Sc);

% The quality measure if two clusters k and k' are merged is given by :
% Q(k,k') = Qinit+2/mn(Sc(k,k')-Card(k,k')). So, we create the matrix
% G = 2/mn(Sc(k,k')-Card(k,k')) and while there are a non diagonal entry
% which is positive, we actually merge the clusters corresponding to its
% indices
G = Sc;

% We want to amalgamate blocks detected during the preprocessing stage and
% those detected during the spectral stage. We do not want to amalgamate
% preprocessing blocks together nor spectral blocks.
assign_blocks = zeros(1,Kc);
assign_blocks(1:nb_bloc_preproc) = 1;
assign_blocks = (assign_blocks==1);
G(assign_blocks,assign_blocks) = -inf;
G(~assign_blocks,~assign_blocks) = -inf;
G(la_diag) = 0;
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
        % so we stop the process.
    elseif(i==j)
        Qref=-1; i=-1;j=-1;
    end
    % If none of the pair can improve the quality measure, we leave the
    % loop and so the function
    if(Qref <= 0)
        Continuer = false;
        % But if we have a couple which improves the quality measure, we
        % amalgame the two initial clusters and we update the variables
    else
        % We have find a couple (i,j) to amalgame, so we put the elements of
        % the jth cluster in the ith cluster, and we remove the jth cluster.
        Kc = Kc-1;
       
        la_diag = 0:Kc:Kc*(Kc-1);
        la_diag = la_diag +[1:Kc];
        % Here is the new quality measure (after amalgamation)
        Qmes = Qmes + 2*Qref;
        % We update the condensed matrix : we add jth column and row to the
        % ith ones and remove them.
        Sc(i,:) = Sc(i,:)+Sc(j,:);
        Sc(:,i) = Sc(i,:)';
        
        Sc(j,:)=[];Sc(:,j) = [];
        G = Sc;
        
        assign_blocks(i)=false;assign_blocks(j) = [];
        G(assign_blocks,assign_blocks) = -inf;
        G(~assign_blocks,~assign_blocks) = -inf;
        G(la_diag) = 0;

        
        % we move the indices of the elements of the new cluster in a same block
        % this can be done since we always have i<j 
        p = [p(1:r(i+1)-1), p(r(j):r(j+1)-1), p(r(i+1):r(j)-1), p(r(j+1):r(end)-1)];
        % we finally update the r vector, to remove jth cluster, and change
        % the number of elements of ith cluster.
        nbelts = r(j+1)-r(j); r(j) = []; r(i+1:j-1) = r(i+1:j-1)+nbelts;
        
    end
end

p = p(:); r = r(:);

end
