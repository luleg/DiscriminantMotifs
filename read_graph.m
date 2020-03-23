path_from_nwk = 'Processing/Networks/'; % path towards folder that contains
% the networks, whose format of networks must be:
% # some infos about the network. Contains as many line as one wants, but
% they have to start with '#' to express this is a comment line
% !n:number of nodes
% !m:numer of edges
% v1 v2 (to indicate that ther is an edge from v1 to v2)
path_from_motifs = 'Processing/CountMotifs/'; % path towards folder that
% contains the decomposition of the networks onto 3-node and 4-node motifs.
% Same format as the motifs function from SNAP

path_to = 'MatNetworks/'; % path toward the folder where the vectors will
% be saved

ext = 'elec';
save_pbm(ext,path_from_nwk,path_from_motifs,path_to)

ext = 'fw';
save_pbm(ext,path_from_nwk,path_from_motifs,path_to)

ext = 'soc';
save_pbm(ext,path_from_nwk,path_from_motifs,path_to)


ext = 'stac';
save_pbm(ext,path_from_nwk,path_from_motifs,path_to)


function [] = save_pbm(ext,path_from_nwk,path_from_motifs,path_to)
% This function is to read all the networks from one field (specified by an
% extension, for instance 'soc' for social networks) that are listed in a
% txt file named liste_ext.txt, and to generate and save a mat structure
% where the networks are loaded (one structure per network) that contains:
% -some informations about the networks, if provided in the initial file,
%   in the field 'entete' (French for header)
% -number of edges (field nb_edges), number of nodes (field nb_nodes)
% -edges (a nb_edgesx2 matrix: if there is an edge v1->v2, then there is a
%   row i in edge such that edges(i,1)=v1, edges(i,2)=v2). Node ids are
%   successive and start with 0.
% -fields motif3, motif4: decomposition of the network onto 3-node and
%   4-node motifs. motif3 is a 13x2 matrix such that
%   motif3(i,:) = [id of ith 3-node motif, nb of this motif within network]
% Parameters :
% ext : the prefixe of the graph
% path_from_network : path towards folder that contains the networks
% path_from_motifs : path towards folder that contains the decomposition of
%                    the networks onto 3-node and 4-node motifs.
% path_to: path toward the folder where the vectors will be saved


nb_reseaux = 0;
lreseaux = {};
fid = fopen(['liste_',ext,'.txt']); % each line in liste_ext.txt is a
% network from field ext
tline = fgetl(fid);
while ischar(tline)
    nb_reseaux = nb_reseaux +1;
    lreseaux{nb_reseaux} = tline; % in lreseaux there is the name of all
    % networks from the field
    tline = fgetl(fid);
end
fclose(fid);

save(['names_',ext],'lreseaux') % these names are loaded in a structure
% saved in a mat file, to be accessed easily

for k = 1:nb_reseaux
    % save_one_file : creates the structure characterising one network
    Pbm = save_one_file(lreseaux{k},path_from_nwk,'.nwk',path_from_motifs,...
        '_motif3-counts.tab','_motif4-counts.tab',1);
    % the structure is saved in a mat file with the same name as the
    % network
    save([path_to,lreseaux{k}],'Pbm')
end


end

function [Pbm] = save_one_file(nwk_name,path_nwk,ext_nwk,path_motif,ext_motif3,ext_motif4,verbose)

if nargin <6
    verbose = 1;
end
if verbose
    disp(nwk_name)
end
fid = fopen([path_nwk,nwk_name,ext_nwk]); % % the networks, whose format of networks must be:
% # some infos about the network. Contains as many line as one wants, but
% they have to start with '#' to express this is a comment line
% !n:number of nodes
% !m:numer of edges
% v1 v2 (to indicate that ther is an edge from v1 to v2)
entete = []; % header in french, keeps information about the network 
% provided in the initial file
tline = fgetl(fid);
while ischar(tline)
    % strings that give some explanations about the network
    while tline(1) == '#' || tline(1) == 'h'
        entete = [entete,tline(2:end),'\n'];
        tline = fgetl(fid);
    end
    % number of nodes and of edges
    while tline(1) == '!'
        if tline(2) == 'n'
            n = str2num(tline(4:end));
            
        else
            m = str2num(tline(4:end));
            edges = zeros(m,2);
            num_edge = 1;
        end
        tline = fgetl(fid);
    end
    edge = split(tline);
    edges(num_edge,1) = str2num(edge{1});
    edges(num_edge,2) = str2num(edge{2});
    num_edge = num_edge + 1;
    if verbose && mod(num_edge,5000) == 0
        disp([nwk_name, ' -- edge number : ', int2str(num_edge)])
    end
    tline = fgetl(fid);
end
fclose(fid);
%%
if verbose
    disp('Processing 3-motif')
end
% motif3, motif4 :
% decomposition of the networks onto 3-node and 4-node motifs.
% Same format as the motifs function from SNAP
motif3 = tdfread([path_motif,nwk_name,ext_motif3]); 
motif3 = [motif3.MotifId, motif3.Count]; % motif3(i,:) = [id of ith 3-node motif, nb of this motif within network]
[~,ind_mot3] = sort(motif3(:,1)); 
motif3 = motif3(ind_mot3,:); % sort to ensure all the networks have the 
% motifs in same order.


if verbose
    disp('Processing 4-motif ')
end

motif4 = tdfread([path_motif,nwk_name,ext_motif4]);
motif4 = [motif4.MotifId, motif4.Count];
[~,ind_mot4] = sort(motif4(:,1));
motif4 = motif4(ind_mot4,:);

% Creation of the matlab structure.
Pbm.entete = entete;
Pbm.nb_nodes = n;
Pbm.nb_edges = m;
Pbm.motif3 = motif3;
Pbm.motif4 = motif4;
Pbm.edges = edges;

end