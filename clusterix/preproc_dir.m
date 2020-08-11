function [Result,nb_block] = preproc_dir(A,min_size)
% PREPROCDIR is the preprocessing stage for CLUSTERIX applying on directed
% graphs. It first ensures the adjacency matrix A has a zero free diagonal by
% adding 1e-8 on all diagonal entries, it then finds a first clustering with
% dmperm and then works on each bi-irreducible blocks separately (these
% blocks coincide with strongly connected component of the graph).
% For each block that is large enough, it scales the matrix thanks with bnwetns
% algorithm (see "A Fast algorithm for matrix balancing" from P.Knight)
% and then seeks the dominant entries (over 0.75) inside these blocks.
% input parameters:
%   A           : the adjacency matrix of a directed graph
%   min_size    : the minimal size for a block to be analysed
% output parameters:
%   Results     : a cell-array of size nb_block that contains useful
% informations about bi-irreducible blocks:
%       -Results{k}.indices : the row and column indices from A that belong to
%       its k-th bi-irreducible block.
%       -If the k-th block B = A(indices,indices) is large enough,
%	T = Results{k}.block is the doubly-stochastic scaling of B.
%       p = Results{k}.perm_row and r = Results{k}.sep_row are such that if
%	length(r)>2, then each row of T(p(r(i):r(i+1)-1,:) contains an entries
%	over 0.75-which will lead to a diag entry over 0.55 in T*transpose(T)-
%	for all i < length(r)-2
%	q = Results{k}.perm_col and s = Results{k}.sep_col play same role for
% columns.
%   nb_block    : the number of bi-irreducible blocks in As


beta = 0.75;
n = size(A,1);
A = sparse(A);
A = A+1e-8*speye(n); % To enforce A to have total support

% If there is more than one bi-irreducible blocks in A, it means that there
% are several strongly connected components in the graph. We will treat these
% components separately (not a problem cause elements that are in different
% components have no chance to finally belong to the same cluster.

[pdmperm,~,rdmperm,~] = dmperm(A); % provides symmetric perm bcs A diag is zero fre
Result = {};
nb_block = length(rdmperm)-1;

for num_block=1:nb_block
    % The indices of the whole matrix that compose the current block
    ind = pdmperm(rdmperm(num_block):rdmperm(num_block+1)-1);
    nk = length(ind);
    % We load them in a cell
    Result{num_block}.indices = ind;
    % We only follow the preprocessing if the block is large enough
    if(length(ind)>=min_size)
        T = A(ind,ind);
        % Scaling of the current bi-irreducible block
        [r,c] = bnewtns(T);
        TT = spdiags(r,0,nk,nk)*T*spdiags(c,0,nk,nk);
        
        % detection of dominant entries
        [er,ec] = find(TT>beta);
        er = unique(er(:)');ec = unique(ec(:)');
        i1 = setdiff(1:nk,er);p = [er,i1];
        j1 = setdiff(1:nk,ec);q = [ec,j1];
        if(isempty(er))
            r = [1 nk+1];
        else
            r = [1:length(er)+1 nk+1];
        end
        if(isempty(ec))
            s = [1 nk+1];
        else
            s = [1:length(ec)+1 nk+1];
        end
        Result{num_block}.block = TT;
        Result{num_block}.perm_row = p;
        Result{num_block}.sep_row  = r;
        Result{num_block}.perm_col = p;
        Result{num_block}.sep_col  = r;
        
        
    end
end

end
