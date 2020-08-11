function [Conf] = compute_scores(classe,fw,elec,stac,soc)
% Compute the precision, recall, f1-score and true positive, false positive
% and false negative elements, as well as the confusion matrix, of the
% classification given in parameter.
% classe(i) = k indicates that the ith network in the test database has
% been labelled as belonging in class k
% fw : vector of indices of the networks in the test database that are
% really foodwebs (same for elec, stac and soc)
% /!\ this code implicitely chooses that label fr fw is 1, label for elec
% is 2, label for stac is 3, and label for soc is 4.
% Return :
% infos : scores to evaluate the classifier quality
% info(1,k) : precision of kth (1 : fw, 2: elec,...) field of networks
% info(2,k) : recall of kth field of networks
% info(3,k) : f1-score for kth field of networks
% info(4,k) : true positive of kth field of networks
% info(5,k) : fasle positive of kth field of networks
% info(6,k) : false negative of kth field of networks
% Conf : confusion matrix :  Conf(i,j) : number of elements from class i 
% that have been labelled as belonging in class j

assign_fw = 1;assign_elec = 2;assign_stac =3;assign_soc = 4;
fw = fw(:)';elec = elec(:)';stac = stac(:)';soc = soc(:)';

Conf = zeros(4);

cls = classe(fw);
[N,edges] = histcounts(cls);
edges = ceil(edges);
Conf(assign_fw,edges(1:end-1)) = N;

cls = classe(elec);
[N,edges] = histcounts(cls);
edges = ceil(edges);
Conf(assign_elec,edges(1:end-1)) = N;

cls = classe(stac);
[N,edges] = histcounts(cls);
edges = ceil(edges);
Conf(assign_stac,edges(1:end-1)) = N;

cls = classe(soc);
[N,edges] = histcounts(cls);
edges = ceil(edges);
Conf(assign_soc,edges(1:end-1)) = N;

end


