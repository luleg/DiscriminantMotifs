function [class_t] = fct_dist2mean(W_a,W_t,class_a)
% This function attributes a class to each indiv in W_t according to its k 
% nearest neighbours
% W_a/t ; [w1,...,wn], where wk is an individual
% W_a is the learning basis, W_t is the test basis:
nbt_reseaux = size(W_t,2);

P = sqdist(W_a,W_t); % P : P(i,j) is the square euclid. distance between 
% the ith indiv from W_a and the jth from W_t

class_t = zeros(1,nbt_reseaux);
for i =1:nbt_reseaux
    col = P(:,i); % column that contains the distance between the ith 
    % network from the test basis and all the networks from the learning 
    % basis
%     [~,ind] = sort(col); 
%     ind = ind(1); % we sort the distances in ascending order and keep 
    %only the nb_k smallest
    [~,ind] = min(col);
    class_t(i)=class_a(ind);
end
end

