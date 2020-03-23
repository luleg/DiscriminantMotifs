%% Motifs of interest
Dt = D(1:nb_ax);
Dt = Dt/sum(Dt); % percentage will be relative to the kept amount of trace

% Construction of Gamma score
Ut = U(:,1:nb_ax); % t first principal axes (Ut : matrix p x nb_ax)
Nu = abs(Ut)'*ones(t_vect,1); % norm1 of each t principal axes
invNu = 1./Nu;
Ut = Ut*diag(invNu); % each column of Ut is now a principal axis with norm1
% equal to 1
Dut = abs(Ut) * diag(Dt); % in Dut a pxt matrix, a coordinate (i,j) is the 
% expression level of parameter i in the principal axis j, relatively to 
% the importance of the axis j (given by Dt(j))
Gamt = Dut*ones(nb_ax,1); % Gamt(i) : Gamma score of the ith canonical
% axis (and hence motif)

if disp_fig
    disp_gamma(Gamt,51,nb_ax,true);
end
