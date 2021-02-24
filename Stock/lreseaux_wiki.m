fid = fopen('liste_wiki_guerre-edition.txt');
lreseaux = {};
cpt=1;
line=fgetl(fid);
while (line~=-1)
    lreseaux{cpt}=line;
    cpt=cpt+1;
    line=fgetl(fid);
end
fclose(fid);
save('names_wiki_guerre-edition','lreseaux')

fid = fopen('liste_wiki_non-neutre.txt');
lreseaux = {};
cpt=1;
line=fgetl(fid);
while (line~=-1)
    lreseaux{cpt}=line;
    cpt=cpt+1;
    line=fgetl(fid);
end
fclose(fid);
save('names_wiki_non-neutre','lreseaux')

fid = fopen('liste_wiki_qualite.txt');
lreseaux = {};
cpt=1;
line=fgetl(fid);
while (line~=-1)
    lreseaux{cpt}=line;
    cpt=cpt+1;
    line=fgetl(fid);
end
fclose(fid);
save('names_wiki_qualite','lreseaux')