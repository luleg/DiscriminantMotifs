fid = fopen('sorted_names.txt'); name_files = {};
tline = fgetl(fid);
i=0;

while(tline~=-1)
    i=i+1;
    name_files{i}=tline;
    tline =fgetl(fid);
end
fclose(fid);

mkdir EmbMat
listing = dir('Embeddings');
nb_files = size(listing,1);
for i=1:nb_files
    
    file = listing(i).name;
    if (length(file)>4 && all(file(end-3:end)=='.csv'))
        underscores = strfind(file,'_');
        end_emb = strfind(file,'.');
        t_emb = file(underscores(end)+1:end_emb-1);
        t_wl = file(underscores(1)+1:underscores(2)-1);
        mkdir(['EmbMat/WL',t_wl,'_Emb',t_emb]);
        
        Emb = csvread(['Embeddings/',file]);
        
        for k = 1:size(Emb,1)
            emb = Emb(k,:)';
            save(['EmbMat/WL',t_wl,'_Emb',t_emb,'/',name_files{k}],'emb');
        end
    end
end
