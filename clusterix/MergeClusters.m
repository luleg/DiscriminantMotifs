function [q, r] = MergeClusters(q1, r1, q2, r2, q0, r0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MERGECLUSTERS overlaps the clustering given by (q1,r1) and
% (q2,r2).
% If (q0,r0) are provided, it also takes apart the corresponding
% blocks and put it as first indices in the returned permutation
% (to always ensure the "mini blocks" from preprocessing are along 
% the top diagonal).


m = length(q1);  G = zeros(m,1);

for k = 1:length(r1)-1
    qk = q1(r1(k):r1(k+1)-1);
    G(qk) = k*m^2*ones(length(qk),1);
end

if (nargin >= 4)
    for k = 1:length(r2)-1
        qk = q2(r2(k):r2(k+1)-1);
        G(qk) = G(qk)+k*m*ones(length(qk),1);
    end
end

if(nargin == 6)
    for k=1:length(r0)-2
        qk = q0(r0(k):r0(k+1)-1);
        G(qk) = (k-1)*ones(length(qk),1);
    end
end

[G,q] = sort(G);
[~,r,~] = unique(G,'first');
r = [r;m+1];

end
