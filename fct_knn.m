function [class_t] = fct_knn(W_a,W_t,class_a)
% This function attributes a class to each indiv in W_t according to its k 
% nearest neighbours
% W_a/t ; [w1,...,wn], where wk is an individual
% W_a is the learning basis, W_t is the test basis:
% class_a :
global nb_k
nbt_reseaux = size(W_t,2);

P = sqdist(W_a,W_t); % P : P(i,j) is the square euclid. distance between 
% the ith indiv from W_a and the jth from W_t

class_t = zeros(1,nbt_reseaux);
for i =1:nbt_reseaux
    col = P(:,i); % column that contains the distance between the ith 
    % network from the test basis and all the networks from the learning 
    % basis
    [~,ind] = sort(col); 
    ind = ind(1:nb_k); % we sort the distances in ascending order and keep 
    %only the nb_k smallest
    class_i = class_a(ind);
    class_i = mode(class_i);
    class_t(i)=class_i;
end
end

function d=sqdist(a,b)
% SQDIST - computes squared Euclidean distance matrix
%          computes a rectangular matrix of pairwise distances
% between points in A (given in columns) and points in B

% NB: very fast implementation taken from Roland Bunschoten

aa = sum(a.*a,1); bb = sum(b.*b,1); ab = a'*b; 
d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end
