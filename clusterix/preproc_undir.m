function [As,Result,nb_block] = preproc_undir(A,min_size)
% PREPROCUNDIR is the preprocessing stage for CLUSTERIX applying on undirected
% graphs. It first ensureis the adjacency matrix A has a zero free diagonal by 
% adding 1e-8 on all diagonal entries, it then scales the matrix thanks to 
% bnwetns algorithm (see "A Fast algorithm for matrix balancing" from P.Knight)
% It finds a first clustering by spearating the bi-irreducible blocks (that 
% coincide with the connected components), with dmperm and then seeks the 
% dominant entries (over 0.55) inside these blocks. 
% input parameters:
%   A 		: the adjacency matrix of an undirected graph
%   min_size	: the minimal size for a block to be analysed
% output parameters:
%   As		: the doubly-stochastic scaling of A+1e-8*Identity
%   Results     : A cell-array of size nb_block that contains useful 
% informations about bi-irreducible blocks:
%	-Results{k}.indices : the row and column indices from As that belong to 
% 	its k-th bi-irreducible block.
% 	-If the k-th block B = As(indices,indices) is large enough,
%	p = Results{k}.perm and r = Results{k}.sep are such that if length(r)>2,
%	then B(p(r(i):r(i+1)-1,p(r(i),r(i+1)-1)  are scalars over 0.55 or matrix 
%       2x2 that have nondiagonal entries over 0.55, for all i < length(r)-2.
%   nb_block 	: the number of bi-irreducible blocks in As 



beta = 0.55;


n = size(A,1);
% Uncomment if you want to reproduce the results of USd
A = sparse(A);
A = A+1e-8*speye(n,n);
% Uncomment if you want to reproduce the results of USw
% A = A+0.15*ones(n);

[r,c] = bnewtns(A);
As = spdiags(r,0,n,n)*A*spdiags(c,0,n,n);

% If there is more than one bi-irreducible blocks in A, it means that there
% are several connected components in the graph. We will treat these
% components separately (this is not a problem because elements that are in different
% components have no chance to finally belong to the same communities).
[p0,~,r0,~] = dmperm(As); % dmperm provides symmetric perm. bcs A diag is zerofree
Result = {};


nb_block = length(r0)-1;

for num_block=1:nb_block
    % The indices of the whole matrix that compose the current block
    ind = p0(r0(num_block):r0(num_block+1)-1);
    nk = length(ind);
    % We load them in a cell
    Result{num_block}.indices = ind;
    % We only follow the preprocessing if the block is large enough
    if(length(ind)>=min_size)
        % Here we take apart dominant elements
        TT = As(ind,ind);
        [er,ec] = find(triu(TT,1)>beta);
        er = er(:)';ec = ec(:)';
        ex = [er;ec]; ex = ex(:)';
        ed = find(diag(TT)>beta);
        ed = ed(:)';
        i1 = setdiff(1:nk,[ed,ex]);
        p = [ed,ex,i1];
        r = [1:length(ed), length(ed)+1:2:length(ed)+1+length(ex)];
        r = [r,r(end)+length(i1)]; r = unique(r);
        Result{num_block}.perm = p;
        Result{num_block}.sep  = r;
    end
end
end

