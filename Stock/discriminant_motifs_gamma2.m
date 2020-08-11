%% Motifs of interest
Dt = D(1:nb_ax);
Dt = Dt/sum(Dt); % percentage will be relative to the kept amount of trace

% Construction of Gamma score
Ut = U(:,1:nb_ax); % t first principal axes (Ut : matrix p x nb_ax)
sq_Ut = Ut.*Ut;

Dut = sq_Ut * diag(Dt); % in Dut a pxt matrix, a coordinate (i,j) is the 
% expression level of parameter i in the principal axis j, relatively to 
% the importance of the axis j (given by Dt(j))
Gamt = Dut*ones(nb_ax,1); % Gamt(i) : Gamma score of the ith canonical
% axis (and hence motif)
