function V = generate_basis(lreseaux,mot3,mot4,mot5)
% function that loads and pre-processes the networks whose names are in the
% matlab structure lreseaux
% return V = [v1,...,vn], where vk is the kth network characterised by its
% decomposition onto 3-and/or 4-and/or 5-node and graphlets

global type_norm
global norm_2

if nargin<2
    mot3=1;
    mot4=1;
    mot5=0;
    t_vect = 199+13;
else
    if ~mot3 && ~mot5
        t_vect = 199;
    elseif ~mot4 && ~mot5
        t_vect = 13;
    elseif ~mot3 && ~mot4
        t_vect = 9364;
    elseif ~mot5
        t_vect = 13+199;
    else
        t_vect = 13+199+9364;
    end
end
where_to = 'Stock/MatNetworks/';

nb_reseaux = length(lreseaux);
V = zeros(t_vect,nb_reseaux);

for i = 1:nb_reseaux
    load([where_to,lreseaux{i}])
    %     disp([num2str(i),'--', lreseaux{i}])
    
    v4 = Pbm.motif4;
    v4 = v4(:,2); %v4 = v4(keep4);
    v3 = Pbm.motif3;
    v3 = v3(:,2); %v3 = v3(keep3);
%     v3(1,:)=;
    if mot5
        v5 = Pbm.motif5;
        v5 = v5(:,2);
    else
        v5 =0;
    end
    
    m = Pbm.nb_edges;
    n = Pbm.nb_nodes;
    if type_norm == 0
        u3 =13/(n*(n-1)*(n-2))*v3;
        u4 =199/(n*(n-1)*(n-2)*(n-3))*v4;
        u5 = 9364/(n*(n-1)*(n-2)*(n-4)*(n-5))*v5;
        
    elseif type_norm == 1
        u3 =13*6/(n*(n-1)*(n-2))*v3;
        u4 =199*24/(n*(n-1)*(n-2)*(n-3))*v4;
        u5 =9364*120/(n*(n-1)*(n-2)*(n-4)*(n-5))*v5;
    elseif type_norm == 2
        u3 =6/(n*(n-1)*(n-2))*v3;
        u4 =24/(n*(n-1)*(n-2)*(n-3))*v4;
        u5 =120/(n*(n-1)*(n-2)*(n-4)*(n-5))*v5;
    elseif type_norm == 3
        u3 =log(13)/(n*(n-1)*(n-2))*v3;
        u4 =log(199)/(n*(n-1)*(n-2)*(n-3))*v4;
        u5 =log(9364)/(n*(n-1)*(n-2)*(n-4)*(n-5))*v5;
    elseif type_norm == 4
        u3 =log(13/(n*(n-1)*(n-2)))*v3;
        u4 =log(199/(n*(n-1)*(n-2)*(n-3)))*v4;
        u5 =log(9364/(n*(n-1)*(n-2)*(n-4)*(n-5)))*v5;
     elseif type_norm == 5
        u3 =v3/norm(v3);
        u4 =v4/norm(v4);
        u5 =v5/norm(v5);
     elseif type_norm == 6
        u3 =sqrt(13)*v3/norm(v3);
        u4 =sqrt(199)*v4/norm(v4);
        u5 =sqrt(9364)*v5/norm(v5);
    elseif type_norm == 7
        u3 =v3;
        u4 =v4;
        u5 =v5;
    end
    u = [];
    if mot3
        u=[u;u3];
    end
    if mot4
        u=[u;u4];
    end
    if mot5
        u=[u;u5];
    end
    %
%     u(setdiff(1:t_vect,[2,1,4,72,18,20,3,14,8]))=0;
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
