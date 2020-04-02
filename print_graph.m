function [G] = print_graph(node,motif,num_fig)
% Takes as an input a motif id, plots and returns the corresponding graph
% Inuts : 
%   - node : number of nodes (3 or 4)
%   - motif : motif id as explained in Sec3 of Supplementary Material
%   - (Optional) num_fig : id of the Figure one wants to plot the graph
% /!\ : It does not check if the id is a valid id.

id = motif;
if node == 3
    power_2 = [0 128 64 32 0 8 4 2 1];
elseif node == 4
    power_2 = [0 16384 8192 4096 2048 1024 512 256 128 64 0 16 8 4 2 1];
end

A = zeros(1,node*node);
cpt = 1;
while motif>0
    tmp = mod(motif,power_2(1));
    power_2(1) = [];
    if tmp < motif
           A(cpt) = 1;
           motif = tmp;
    end
    cpt = cpt+1;
end

A = fliplr(A);
A = reshape(A,[node node]);
A = A';
G = digraph(A);
if nargin <3, num_fig = 1000;end
figure(num_fig),clf
pl = plot(G);


if node == 3
    pl.XData = [0 -1 1];
    pl.YData = [1 0 0];
elseif node == 4
    pl.XData = [0 1 1 0];
    pl.YData = [1 1 0 0];
end
ttl = title(['id : ',int2str(id)]);
ttl.FontSize = 20;
pl.MarkerSize = 20;
pl.NodeFontSize = 30;
pl.NodeFontWeight = 'bold';

pl.LineWidth = 5;
pl.ArrowSize = 30;

end

