function [p] = clusterix(A,nsteps,deep,min_taille,is_directed)
% CLUSTERIX finds the underlying block structure of a square matrix A essentially
% with spectral partitioning. It is a tool used to cluster datasets, then A is
% expected to be the adjacency matrix of a directed or undirected graph.
% input parameters:
%   A           : the matrix to cluster
%   nsteps      : maximum number of algorithm global iterations
%   deep        : number of singular vectors analysed by global iteration
%   min_taille  : minimum size of bi-irreducible block to analyse
%   is_directed : true if A is the adjancecy matrix of a directed graph, false
% otherwise (generally if A is symmetric)
% output parameter :
%   p           : assignation of an entry to a block (p(i) = j means that the
% i-th node of the dataset) belongs to the j-th cluster
n = size(A,1);

% Set up of parameters
switch nargin
    case 1
        nsteps = 3;
        deep = 2;
        min_taille = floor(n/10);
        is_directed = issymmetric(A);
    case 2
        deep = 2;
        min_taille = floor(n/10);
        is_directed = issymmetric(A);
    case 3
        min_taille = floor(n/10);
        is_directed = issymmetric(A);
    case 4
        is_directed = issymmetric(A);
end

if ~is_directed && ~issymmetric(A)
    A = 1/2*(A+A');
end

% Preprocessing : finds the bi-irreducible blocks in A, scales them, separates
% dominant entries in the scaled blocks

if ~is_directed
    [T,Result,nb_block] = preproc_undir(A,min_taille);
else
    [Result,nb_block] = preproc_dir(A,min_taille);
end


pp = [];
rr = 1;


for k = 1:nb_block
    ind = Result{k}.indices;ind = ind(:);
    nk = length(ind);
    if nk >=min_taille % We only analyse a block if its size is large enough
        if is_directed
            As = Result{k}.block;
            p = Result{k}.perm_row; r = Result{k}.sep_row; p = p(:); r = r(:);
            q = Result{k}.perm_col; s = Result{k}.sep_col; q = q(:); s = s(:);
            nb_block_row = length(r)-2; nb_block_col = length(s)-2;
            [p,r] = CoarseClusterix(As*As',p,r,deep,nsteps); % the main process
            [q,s] = CoarseClusterix(As'*As,q,s,deep,nsteps); % the main process
            [p,r] = PostProcDir(As,p,r,nb_block_row,q,s,nb_block_col);
            p = p(:); r = r(:); q = q(:); s = s(:);
        else
            As = T(ind,ind);
       
            p = Result{k}.perm; p = p(:); r = Result{k}.sep;r = r(:);
            nb_bloc_preproc = length(r)-2;
            [p,r] = CoarseClusterix(As,p,r,deep,nsteps); % the main process
            [p,r] = PostProcUnDir(As,p,r,nb_bloc_preproc);
            p = p(:);r = r(:);
        end
    else
        p = (1:nk)';
        r = [1 nk+1]';
    end
    pp = [pp;ind(p)];
    rr = [rr;r(2:end)+rr(end)-1];
end

% creation of the assignation vector from pp and rr
p = zeros(1,n);
for i = 1:length(rr)-1
    p(pp(rr(i):rr(i+1)-1)) = i;
end

end
