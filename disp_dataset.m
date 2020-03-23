function [] = disp_dataset(W,fw,elec,stac,soc,fig)
% Display the individuals projected onto the 3 main axes, and onto each
% unique component kept to obatin the target percentage of trace
threshold = size(W,2);
nb_fw = length(fw);nb_elec = length(elec);
nb_stac = length(stac);nb_soc = length(soc);
figure(fig), clf, hold on, grid on
title('Projection onto the 3 principal axes')
if (nb_fw>0), plot3(W(fw,1),W(fw,2),W(fw,3),'b+','markersize',7,'linewidth',2);end
if (nb_elec>0), plot3(W(elec,1),W(elec,2),W(elec,3),'mo','markersize',7,'linewidth',2);end
if (nb_stac>0), plot3(W(stac,1),W(stac,2),W(stac,3),'kx','markersize',7,'linewidth',2);end
if (nb_soc>0), plot3(W(soc,1),W(soc,2),W(soc,3),'r*','markersize',7,'linewidth',2);end

lgd = legend('foodweb','elec. circ','disc. struc','social');

for ax = 1:threshold
    figure(10*fig+ax-1),clf,hold on
    title(['Projection onto the ',int2str(ax),'th axis'])
    if (nb_fw>0),plot(W(fw,ax),zeros(nb_fw,1),'b+','markersize',7,'linewidth',1);end
    if (nb_elec>0),plot(W(elec,ax),zeros(nb_elec,1),'mo','markersize',7,'linewidth',1);end
    if (nb_stac>0),plot(W(stac,ax),zeros(nb_stac,1),'kx','markersize',7,'linewidth',1);end
    if (nb_soc>0),plot(W(soc,ax),zeros(nb_soc,1),'r*','markersize',7,'linewidth',1);end
end


end

