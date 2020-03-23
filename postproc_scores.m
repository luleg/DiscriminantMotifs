function [tab_prec,tab_rec,tab_f1,Conf] = postproc_socres(tens_infos,tens_Conf)
% Postprocessing of the quality scores: from tens_infos and tens_Conf,
% return 3 tables that contains the minimum, maximum, mean and standard 
% deviation over all the iterations of respectively precision, recall, and
% f1-score, for each field independantly. It also compute the mean
% normalised Confusion matrix: in Conf(i,j) is the rate of elements from 
% class i that have been classified in class j, in average, overall the 
% iterationas .


global nb_ites

prec = zeros(nb_ites,4);
recall = zeros(nb_ites,4);
f1 = zeros(nb_ites,4);

Conf = tens_Conf(:,:,1);
D = Conf*ones(4,1);
D = diag(1./D);
Conf = zeros(4);
for i =1:nb_ites
    Conf = Conf+D*tens_Conf(:,:,i);
    prec(i,:) = tens_infos(1,:,i);
    recall(i,:) = tens_infos(2,:,i);
    f1(i,:) = tens_infos(3,:,i);
end

Conf = 1/nb_ites*Conf;

tab_prec = [min(prec(:,1)), min(prec(:,2)), min(prec(:,3)), min(prec(:,4));
    max(prec(:,1)), max(prec(:,2)), max(prec(:,3)), max(prec(:,4));
    mean(prec,1);
    std(prec,0,1)];

tab_rec = [min(recall(:,1)), min(recall(:,2)), min(recall(:,3)), min(recall(:,4));
    max(recall(:,1)), max(recall(:,2)), max(recall(:,3)), max(recall(:,4));
    mean(recall,1);
    std(recall,0,1)];

tab_f1 = [min(f1(:,1)), min(f1(:,2)), min(f1(:,3)), min(f1(:,4));
    max(f1(:,1)), max(f1(:,2)), max(f1(:,3)), max(f1(:,4));
    mean(f1,1);
    std(f1,0,1)];




end

