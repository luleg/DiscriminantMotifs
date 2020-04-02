% To generate the Table of Sec 2.2 that contains statistics about networks
clear

% First food webs are processed
load('names_fw')
nb_reseaux = length(lreseaux);

stats_fw = zeros(3,nb_reseaux); % this table will contained info about each
% food web network

for i = 1:nb_reseaux % For each food web in the benchmark
    load(['MatNetworks/',lreseaux{i}])
    disp([num2str(i),'--', lreseaux{i}])
    m = Pbm.nb_edges; % its number of edges
    n = Pbm.nb_nodes; % its number of nodes
    stats_fw(:,i) = [n;m;m/n]; % saving n,m,m/n in the corresponding column
    % of stats_fw
end

stats_fw1 = zeros(3,5); % will contains stats about the food webs dataset:
% namely:
% 1) max value,
% 2) min value
% 3) mean value
% 4) median
% 5) standard deviation

% of the number of nodes on 1st row
stats_fw1(1,:) = [max(stats_fw(1,:)), min(stats_fw(1,:)), mean(stats_fw(1,:)), median(stats_fw(1,:)), std(stats_fw(1,:))];
% of the number of edges on 2nd row
stats_fw1(2,:) = [max(stats_fw(2,:)), min(stats_fw(2,:)), mean(stats_fw(2,:)), median(stats_fw(2,:)), std(stats_fw(2,:))];
% of the density on 3rd row
stats_fw1(3,:) = [max(stats_fw(3,:)), min(stats_fw(3,:)), mean(stats_fw(3,:)), median(stats_fw(3,:)), std(stats_fw(3,:))];

% same as food webs for electronic circuits :
load('names_elec')
nb_reseaux = length(lreseaux);

stats_elec = zeros(3,nb_reseaux);

for i = 1:nb_reseaux
    load(['MatNetworks/',lreseaux{i}])
    disp([num2str(i),'--', lreseaux{i}])
    m = Pbm.nb_edges;
    n = Pbm.nb_nodes;
    stats_elec(:,i) = [n;m;m/n];
end

stats_elec1 = zeros(3,5);
stats_elec1(1,:) = [max(stats_elec(1,:)), min(stats_elec(1,:)), mean(stats_elec(1,:)), median(stats_elec(1,:)), std(stats_elec(1,:))];
stats_elec1(2,:) = [max(stats_elec(2,:)), min(stats_elec(2,:)), mean(stats_elec(2,:)), median(stats_elec(2,:)), std(stats_elec(2,:))];
stats_elec1(3,:) = [max(stats_elec(3,:)), min(stats_elec(3,:)), mean(stats_elec(3,:)), median(stats_elec(3,:)), std(stats_elec(3,:))];


% same as food webs for discourse structure :
load('names_stac')
nb_reseaux = length(lreseaux);

stats_stac = zeros(3,nb_reseaux);

for i = 1:nb_reseaux
    load(['MatNetworks/',lreseaux{i}])
    disp([num2str(i),'--', lreseaux{i}])
    m = Pbm.nb_edges;
    n = Pbm.nb_nodes;
    stats_stac(:,i) = [n;m;m/n];
end

stats_stac1 = zeros(3,5);
stats_stac1(1,:) = [max(stats_stac(1,:)), min(stats_stac(1,:)), mean(stats_stac(1,:)), median(stats_stac(1,:)), std(stats_stac(1,:))];
stats_stac1(2,:) = [max(stats_stac(2,:)), min(stats_stac(2,:)), mean(stats_stac(2,:)), median(stats_stac(2,:)), std(stats_stac(2,:))];
stats_stac1(3,:) = [max(stats_stac(3,:)), min(stats_stac(3,:)), mean(stats_stac(3,:)), median(stats_stac(3,:)), std(stats_stac(3,:))];


% same as food webs for social networks :
load('names_soc')
nb_reseaux = length(lreseaux);

stats_soc = zeros(3,nb_reseaux);

for i = 1:nb_reseaux
    load(['MatNetworks/',lreseaux{i}])
    disp([num2str(i),'--', lreseaux{i}])
    m = Pbm.nb_edges;
    n = Pbm.nb_nodes;
    stats_soc(:,i) = [n;m;m/n];
end

stats_soc1 = zeros(3,5);
stats_soc1(1,:) = [max(stats_soc(1,:)), min(stats_soc(1,:)), mean(stats_soc(1,:)), median(stats_soc(1,:)), std(stats_soc(1,:))];
stats_soc1(2,:) = [max(stats_soc(2,:)), min(stats_soc(2,:)), mean(stats_soc(2,:)), median(stats_soc(2,:)), std(stats_soc(2,:))];
stats_soc1(3,:) = [max(stats_soc(3,:)), min(stats_soc(3,:)), mean(stats_soc(3,:)), median(stats_soc(3,:)), std(stats_soc(3,:))];


% for the whole datase :
stats_all = [stats_fw,stats_elec,stats_stac,stats_soc];

stats_all1 = zeros(3,5);
stats_all1(1,:) = [max(stats_all(1,:)), min(stats_all(1,:)), mean(stats_all(1,:)), median(stats_all(1,:)), std(stats_all(1,:))];
stats_all1(2,:) = [max(stats_all(2,:)), min(stats_all(2,:)), mean(stats_all(2,:)), median(stats_all(2,:)), std(stats_all(2,:))];
stats_all1(3,:) = [max(stats_all(3,:)), min(stats_all(3,:)), mean(stats_all(3,:)), median(stats_all(3,:)), std(stats_all(3,:))];