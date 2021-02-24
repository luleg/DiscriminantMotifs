ext_nwk = 'TriadCencus.txt';
ext_rdm = 'ExpectedTriadFreqWNumNodeEdge.txt';
epsi = 4;
mkdir EmbMat
t = cputime;
load ../Stock/names_fw

for i = 1:length(lreseaux)
    nwk = lreseaux{i};
    fid = fopen(['MotifCounts/',nwk,ext_nwk]);
    counts_nwk = fscanf(fid,'%lf');
    fclose(fid);
    fid = fopen(['MotifCounts/',nwk,ext_rdm]);
    counts_rdm = fscanf(fid,'%lf');
    fclose(fid);
    emb = (counts_nwk-counts_rdm)./(counts_nwk+counts_rdm+epsi);
    emb = emb/norm(emb,2);
    save(['EmbMat/',nwk],'emb');
end


load ../Stock/names_elec

for i = 1:length(lreseaux)
    nwk = lreseaux{i};
    fid = fopen(['MotifCounts/',nwk,ext_nwk]);
    counts_nwk = fscanf(fid,'%lf');
    fclose(fid);
    fid = fopen(['MotifCounts/',nwk,ext_rdm]);
    counts_rdm = fscanf(fid,'%lf');
    fclose(fid);
    emb = (counts_nwk-counts_rdm)./(counts_nwk+counts_rdm+epsi);
    emb = emb/norm(emb,2);
    save(['EmbMat/',nwk],'emb');
end

load ../Stock/names_stac

for i = 1:length(lreseaux)
    nwk = lreseaux{i};
    fid = fopen(['MotifCounts/',nwk,ext_nwk]);
    counts_nwk = fscanf(fid,'%lf');
    fclose(fid);
    fid = fopen(['MotifCounts/',nwk,ext_rdm]);
    counts_rdm = fscanf(fid,'%lf');
    fclose(fid);
    emb = (counts_nwk-counts_rdm)./(counts_nwk+counts_rdm+epsi);
    emb = emb/norm(emb,2);
    save(['EmbMat/',nwk],'emb');
end

load ../Stock/names_soc

for i = 1:length(lreseaux)
    nwk = lreseaux{i};
    fid = fopen(['MotifCounts/',nwk,ext_nwk]);
    counts_nwk = fscanf(fid,'%lf');
    fclose(fid);
    fid = fopen(['MotifCounts/',nwk,ext_rdm]);
    counts_rdm = fscanf(fid,'%lf');
    fclose(fid);
    emb = (counts_nwk-counts_rdm)./(counts_nwk+counts_rdm+epsi);
    emb = emb/norm(emb,2);
    save(['EmbMat/',nwk],'emb');
end
t = cputime-t;