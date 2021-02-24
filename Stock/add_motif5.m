fold_mot5='Processing/Motif5/';
fold_matnwk= 'MatNetworks/';
nwks = dir(fold_mot5);

for num_nwk = 1:size(nwks,1)
    nwk = nwks(num_nwk).name;
    if ~contains(nwk,'.mat')
        continue
    end
    load([fold_mot5,nwk])
    try
        load([fold_matnwk,nwk])
        Pbm.motif5 = motif5;
        save([fold_matnwk,nwk],'Pbm')
    catch
        disp(['Pbm with file ', nwk])
    end
    
end
