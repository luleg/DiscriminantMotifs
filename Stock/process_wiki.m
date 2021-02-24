fold = 'Processing/Wiki_Graph_remove_40_0_500/Nwks/';
foldMot = 'Processing/Wiki_Graph_remove_40_0_500/Motifs/';
files=dir(fold);

for file =1:size(files,1)
    file = files(file).name;
    if ~contains(file,'wiki')
        continue
    end

    E = readmatrix([fold,file]);
    Pbm = struct;
    Pbm.nb_nodes = E(1,1);
    Pbm.nb_edges= E(1,2);
    Pbm.edges=E;
    fmot3 = [file(1:end-4),'_3.txt'];
    fmot4 = [file(1:end-4),'_4.txt'];
    mot3 = readmatrix([foldMot,fmot3]);
    mot4 = readmatrix([foldMot,fmot4]);
    [~,inds] = sort(mot3(:,1));
    mot3 = mot3(inds,:);
    [~,inds] = sort(mot4(:,1));
    mot4 = mot4(inds,:);
    Pbm.motif3 = mot3;
    Pbm.motif4 = mot4;
    save(['MatNetworks/',file(1:end-4)],'Pbm');
end