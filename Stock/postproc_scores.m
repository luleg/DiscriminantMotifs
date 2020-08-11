function [tab_prec,tab_rec,tab_f1,Conf] = postproc_scores(tens_Conf)
% Postprocessing of the quality scores: from tens_infos and tens_Conf,
% return 3 tables that contains the minimum, maximum, mean and standard 
% deviation over all the iterations of respectively precision, recall, and
% f1-score, for each field independantly. It also compute the mean
% normalised Confusion matrix: in Conf(i,j) is the rate of elements from 
% class i that have been classified in class j, in average, overall the 
% iterationas .

nb_ites = size(tens_Conf,3);


prec = zeros(nb_ites,4);
rec = zeros(nb_ites,4);
f1 = zeros(nb_ites,4);
CConf = zeros(4);

e = ones(4,1);

D = tens_Conf(:,:,1)*e;
D = diag(1./D);

for i=1:nb_ites
    Conf = tens_Conf(:,:,i);
    tp = diag(Conf); tp = tp(:);
    fn = Conf*e-tp;
    fp = (e'*Conf)'-tp;
    prec(i,:) = tp./(tp+fp);
    rec(i,:) = tp./(tp+fn);
    f1(i,:) = 2*(prec(i,:).*rec(i,:))./(prec(i,:)+rec(i,:));
    CConf = CConf+D*Conf;
end

Conf = 1/nb_ites*CConf;

tab_prec = [min(prec(:,1)), min(prec(:,2)), min(prec(:,3)), min(prec(:,4));
    max(prec(:,1)), max(prec(:,2)), max(prec(:,3)), max(prec(:,4));
    mean(prec,1);
    std(prec,0,1)];

tab_rec = [min(rec(:,1)), min(rec(:,2)), min(rec(:,3)), min(rec(:,4));
    max(rec(:,1)), max(rec(:,2)), max(rec(:,3)), max(rec(:,4));
    mean(rec,1);
    std(rec,0,1)];

tab_f1 = [min(f1(:,1)), min(f1(:,2)), min(f1(:,3)), min(f1(:,4));
    max(f1(:,1)), max(f1(:,2)), max(f1(:,3)), max(f1(:,4));
    mean(f1,1);
    std(f1,0,1)];




end

