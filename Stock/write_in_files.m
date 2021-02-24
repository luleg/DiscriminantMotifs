load RandPermDatasets.mat
str_elec = "";
str_fw = "";
str_soc = "";
str_stac = "";
nb_train = 40;

for i=1:500
    tmp_elec = string(zeros(1,nb_train));
    load names_elec.mat
    for k = 1:nb_train
        tmp_elec(k) = lreseaux{lind_elec(i,k)};
    end
    str_elec = str_elec+join(tmp_elec,',')+newline;
    
    tmp_fw = string(zeros(1,nb_train));
    load names_fw.mat
    for k = 1:nb_train
        tmp_fw(k) = lreseaux{lind_fw(i,k)};
    end
    str_fw = str_fw+join(tmp_fw,',')+newline;
    
    tmp_soc = string(zeros(1,nb_train));
    load names_soc.mat
    for k = 1:nb_train
        tmp_soc(k) = lreseaux{lind_soc(i,k)};
    end
    str_soc = str_soc+join(tmp_soc,',')+newline;
    
    tmp_stac = string(zeros(1,nb_train));
    load names_stac.mat
    for k = 1:nb_train
        tmp_stac(k) = lreseaux{lind_stac(i,k)};
    end
    str_stac = str_stac+join(tmp_stac,',')+newline;
end

fid = fopen('perms_elec.txt','w');
fprintf(fid,'%s',str_elec)
fclose(fid)
fid = fopen('perms_fw.txt','w');
fprintf(fid,'%s',str_fw)
fclose(fid)
fid = fopen('perms_soc.txt','w');
fprintf(fid,'%s',str_soc)
fclose(fid)
fid = fopen('perms_stac.txt','w');
fprintf(fid,'%s',str_stac)
fclose(fid)