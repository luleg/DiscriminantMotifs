function V = generate_basis(lreseaux)
% function that loads and pre-processes the networks whose names are in the
% matlab structure lreseaux
% return V = [v1,...,vn], where vk is the kth network characterised by its
% decomposition onto 3-node and 4-node motifs
global t_vect
global type_norm
global norm_2

nb_reseaux = length(lreseaux);
V = zeros(t_vect,nb_reseaux);

for i = 1:nb_reseaux
    load(['MatNetworks/',lreseaux{i}])
    %     disp([num2str(i),'--', lreseaux{i}])
    v4 = Pbm.motif4;
    v4 = v4(:,2);
    v3 = Pbm.motif3;
    v3 = v3(:,2);
    m = Pbm.nb_edges;
    n = Pbm.nb_nodes;
    if type_norm == 0
        u = [13*v3/(n*(n-1)*(n-2));199*v4/(n*(n-1)*(n-2)*(n-3))];
    elseif type_norm == 1
        u = [13*6*v3/(n*(n-1)*(n-2));199*24*v4/(n*(n-1)*(n-2)*(n-3))];
    elseif type_norm == 2
        u = [13*v3/norm(v3,1);199*v4/norm(v4,1)];
    end
    %
    if norm_2 == 1
        V(:,i) = u/norm(u);
    else
        V(:,i) = u;
    end
end

if norm_2 == 2
    for i=1:t_vect
        V(i,:) = V(i,:)/norm(V(i,:));
    end
    for i=1:nb_reseaux
        V(:,i) = V(:,i)/norm(V(:,i));
    end
    
end

end
