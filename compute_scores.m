function [Conf,infos] = compute_scores(classe,fw,elec,stac,soc)
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

fw = fw(:)';elec = elec(:)';stac = stac(:)';soc = soc(:)';

Conf = zeros(4);
infos= zeros(6,4);

[row_confu,col_infos]=partial_score(classe,1,4,fw,[elec,stac,soc]);
infos(:,1) = col_infos;
Conf(1,:) = row_confu;

[row_confu,col_infos]=partial_score(classe,2,4,elec,[fw,stac,soc]);
infos(:,2) = col_infos;
Conf(2,:) = row_confu;

[row_confu,col_infos]=partial_score(classe,3,4,stac,[fw,elec,soc]);
infos(:,3) = col_infos;
Conf(3,:) = row_confu;

[row_confu,col_infos]=partial_score(classe,4,4,soc,[fw,elec,stac]);
infos(:,4) = col_infos;
Conf(4,:) = row_confu;
end

function [row_confu,col_infos]=partial_score(classe,lab_class,nb_class,indices_in,indices_out)

row_confu = zeros(1,nb_class);
col_infos = zeros(6,1);

cl = classe(indices_in);
ncl = classe(indices_out);
tp = sum(cl == lab_class);
fn = sum(cl~=lab_class);
fp = sum(ncl==lab_class);
prec = tp/(tp+fp); % Precision
rec = tp/(tp+fn); % Recall
f1 = 2*prec*rec/(prec+rec); % F1-score
col_infos(1)= prec;  
col_infos(2)= rec;
col_infos(3) = f1;
col_infos(4) = tp;
col_infos(5) = fp;
col_infos(6) = fn;

for i=1:length(indices_in)
    
    row_confu(cl(i)) = row_confu(cl(i))+1;
end

end

