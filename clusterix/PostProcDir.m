function [ps,rs] = PostProcDir(A,p,r,nb_bloc_row,q,s,nb_bloc_col)
% POSTPROCDIR is the postprocessing of CLUSTERIX for undirected graph.
% It first applies the postproc for undirected graphs on both normal
% equations.
% Then, it overlaps the clustering for rows and the one for columns and
% applies the amalgamation process on this overlapping.
% INPUT
%    A  : the doubly-stochastic scaling of the adjacency matrix of the
% (sub)graph we want to clusterize.
%    p,r : the row clustering (p(r(i):r(i+1)-1) are the indices of the rows
% that belong to the i-th cluster)
%    nb_bloc_row : the number of row clusters taken apart during the
% preprocessing stage
%    q,s,nb_bloc_col : same thing for the columns
% OUTPUT
%    ps,rs : the final clustering

m = length(p);

AAt = A*A';
AtA = A'*A;
[p,r] = PostProcUnDir(AAt,p,r,nb_bloc_row);
[q,s] = PostProcUnDir(AtA,q,s,nb_bloc_col);


[ps, rs] = MergeClusters(p,r, q,s );

A =1/2*(A+A');
[ps,rs] = Amalgamation(A,ps,rs);

end
