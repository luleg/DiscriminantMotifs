function [p1,r1] = CoarseClusterix(A,p,r,deep,nsteps)
% COARSECLUSTERIX finds the underlying block diagonal structure of
% a symmetric stochastic matrix A by means of spectral clustering.
% input parameters:
%    A	     : the matrix we want to detect the block structure
%    p	     : a permutation of rows and columns of A such that
% dominant entries of A(p,p) are along the top diagonal.
%    r	     : A(p(r(i):r(i+1)-1,p(r(i):r(i+1)-1)) is a scalar over
%  0.55 or a 2x2 matrix with non diagonal entries over 0.55.
%    deep    : number of eigenvectors analysed per iteration
%    nsteps  : max number of iterations


% tol is used as a tolerance to the computation of singular vectors
tol = 1e-8;
% tol_ratio used as a threshold for convergence
tol_ratio = 1e-2;

m = size(A,1);

p0 = p; p0 = p0(:); r0 = r;
p1 = p0;r1 = r0;
Qmes = QualityMeasure( A, p1, r1 );
Qmes_l(1) = Qmes;
nbClust_l(1) = 1;

% Creation of the ONB of the blocks, to project into the orthogonal
% subspace
Kr1 = length(r1)-1;
Vr = zeros(m,Kr1);
for k = 1:Kr1
    Vr(p1(r1(k):r1(k+1)-1), k) = 1/sqrt(r1(k+1)-r1(k));
end

cv = false;
% options for spectral calculations.
opts.tol=tol; opts.maxit=600; opts.disp=0;opts.issym = true;
% For Canny filter width:
maxwidth = (r0(end)-r0(end-1))/10;
for st = 1:nsteps
    if(~cv)
        % Computation of the new set of leading eigenvectors
        defr = @(x) ( A*(x - Vr*(Vr'*x)) - Vr*(Vr'*(A*(x - Vr*(Vr'*x)))) );
        [U,E,flag] = eigs(defr,m,deep,'LA',opts);  if flag, fprintf('eigs failed to converge'), end
        [~,JJ] = sort(diag(E),1,'descend'); U = U(:,JJ);
        Ur = U(:,1:deep)*diag(2*(U(1,1:deep)>=0)-1);

        % Analysis of the steps in this vectors with Canny filter and edge
        % refinement process.
        [pp,rr] = StepDetection(m,Ur,maxwidth);
        % Overlapping with the clustering found during previous iterations
        [pp,rr] = MergeClusters(p1,r1,pp,rr,p0,r0);

        % Amalgamation process: we do not want to amalgamate dominant entries with the rest of the blocks,
	% so we keep them apart
        paux = pp(r0(end-1):r0(end)-1);
        raux = rr-r0(end-1)+1;raux = raux(raux>0);
        Aaux = A(paux,paux);
        [indaux,raux] = Amalgamation( Aaux,1:length(paux),raux);
        pp = [pp(r0(1):r0(end-1)-1); paux(indaux)];
        rr = [r0(1:end-1); raux(2:end)-1+r0(end-1)];
	% Modularity of the clustering at the end of the iteration
        Qmes = QualityMeasure( A, pp, rr );

        % Stopping criterion
        Kr1b = length(rr)-1;
        Vb = zeros(m,Kr1b);
        for k = 1:Kr1b
            Vb(rr(k):rr(k+1)-1, k) = 1/sqrt(rr(k+1)-rr(k));
        end
        Vb(pp,:) = Vb;
        JJ = Vb*Vb';

        if(Kr1b==Kr1 && all(all(sortrows(Vb')==sortrows(Vr'))))
%           disp('cv');
            Qmes_l=Qmes_l(1:st);
            nbClust_l = nbClust_l(1:st);
            cv = true;
        elseif(Qmes<Qmes_l(st))
            pp = p1;rr = r1;
            Qmes_l(st+1) = Qmes_l(st); nbClust_l(st+1) = nbClust_l(st);
 %          disp('cv : Qmes_crt < Qmes_ref')
            Qmes_l=Qmes_l(1:st);
            nbClust_l = nbClust_l(1:st);
            cv = true;
            % if the modu is better for the new clustering, and if the number
            % of clusters is smaller, of course, we accept the new clustering
        elseif(length(rr)-1<nbClust_l(st))
            Qmes_l(st+1) = Qmes; nbClust_l(st+1) = length(rr)-1;
%            disp('Qmes_crt > Qmes_ref et nbClust_crt < nbClust_ref')
            Vr = Vb;
            % if the modu of the new clustering is better, and the number of
            % clusters is greater, we compute the ratio : (real gain)/(ideal case gain)
        else
            ratio = (Qmes-Qmes_l(st))/(1/nbClust_l(st)-1/(length(rr)-1));
            % if the ratio is over the thershold tol_ratio we accept the
            % new clustering
            if(ratio>tol_ratio)
                Qmes_l(st+1) = Qmes;nbClust_l(st+1) = length(rr)-1;
 %               disp('Qmes_crt > Qmes_ref et ratio > tol_ratio')
                Vr = Vb;
                % else we reject it
            else
                Qmes_l(st+1) = Qmes_l(st); nbClust_l(st+1) = nbClust_l(st);
  %              disp('cv : Qmes_crt > Qmes_ref et ratio < tol_ratio')
                Qmes_l=Qmes_l(1:st);
                nbClust_l = nbClust_l(1:st);
                pp = p1;rr = r1;
                cv = true;
            end
        end
        % else we update the parameters
        Kr1 = size(Vr,2);
        p1 = pp;  r1 = rr;
    end

end
end


%% Quality Measure
function Out = QualityMeasure( As, p1, r1 )
gamma = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the quality measure for the given clusterisating.
%
% INPUT :
% As : the matrix to clusterize
% p1 : contains the reorganization of the rows of As to make appear row
%   clusters
% r1 : delimitation of the row cluster : p1(r1(k):r1(k+1)-1) contains the
%   indices of the rows in As which belongs to the kth row cluster
%
% OUTPUT :
% Out : contains the quality measure for the row clustering.

m = size(As,1);
Kc1 = length(r1)-1;
As=As(p1,p1);

U = zeros(m,Kc1);
M = zeros(Kc1,1);
for k = 1:Kc1
    M(k) = r1(k+1)-r1(k);
    U(r1(k):r1(k+1)-1, k) = 1;
end

E = U'*As*U;
Out = sum(diag(E)/m - gamma*((M/m).^2));
end


%% Step Detection
function [p1,r1] = StepDetection(m,U,maxwidth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEPDETECTION analyses in turn the sorted leading eigenvectors of A
% (stored as columns of U) to find their steps, by considering them as a
% signal and applying the Canny filter for different size of width on them.
% It detects the resulting maxima with searchMax function.
% It then submitt this clustering to the EdgeRefinement process to ensure
% the steps are sharp enough, and then merges the resulted clustering with
% those obtained by means of previous vectors.
% input :
%	m : length of the vectors
%	U : each column of U is a vector to analyse
%	maxwidth : the maximal size for Canny filter width
% output :
%	p1, r1 : define the resulting clustering: p(r(i):r(i+1)-1) are a
% subset of 1:m that contains the the indices of one block of A

pref = 1:m;rref=[1 m+1];
m1aux = 0;
deep = size(U,2);

% Centre and normalise the coordinates of the vectors
Norme = diag(U'*U);Norme = sqrt(Norme);Norme=1./Norme; U = U*diag(Norme);

% size of slidding windows for Canny filter
width(1) = max(round(maxwidth),3);
width(2) = max(round(maxwidth/5),3);

for nbvs=1:deep

    u1 = U(:,nbvs);
    % Order vector in ascending order
    [u,p] = sort(u1,'ascend');

    % Step Detection
    %  dDatau = convolution with 1st derivative Gaussian
    %  width = width of slidding window

    dDatauSum = 0;

    % Application of Canny filter on the vector for each window size
    for i=1:length(width)
        dDatau = CreateGaussScaleSpace(u, width(i));
        dDatauSum = dDatauSum + dDatau;

    end
    X = dDatauSum;
    % Maximum detection
    ru = SearchMax(X);

    [p1,r1,m1] = EdgeRefinementV4(U,ru,nbvs);
%     figure(10),clf
%     plot(u1(p1))
%     for kki = 2:length(r1)-1
%         figure(10), hold on,
%         plot([r1(kki)-1/2,r1(kki)-1/2],[min(u1),max(u1)],'k-')
%     end
%     m1
%     pause
    % if we are in the case where m1<0, it means that the edges provided by
    % the call to Affich_sep are really consistant.
    if(m1<0)
        % By the way, if the clusterings found with previous singular vectors
        % are not that consistant, we forgot them.
        if(m1aux>0)
            rref = [1 m+1];
        end
        % and we warn the process that the current clustering is consistant
        m1aux=-1;
        % if we are in the case where the clustering provided by the call to
        % Affich_sep is not that consistant, we only keep the clustering iff
        % the edge inside is sharper than the one in the clustering found
        % before
    elseif((m1>=0 && m1aux<0)||(m1<=m1aux))
        r1 = [1 m+1];
    else
        rref = [1 m+1];
        m1aux = m1;
    end
    % Overlapping witgh the previous found clustering
    [p1, r1] = MergeClusters(p1, r1, pref, rref);
    pref=p1;rref=r1;
end
p1 = pref;r1 =rref;

end


%% SearchMax
function [ru] = SearchMax(dDatau)
% Input :
%  dDatau : convolution product with the 1st Gaussian derivative with the eigenvector

% Output :
% ru  : step interval definition of the eigenvector
tolmax = 1e-8;

m = length(dDatau);
% Local maximum and minimum detection for u
aux =  [dDatau(:)';1:m];
ind_a_sup = find(abs(diff(dDatau))<=tolmax);

aux(:,ind_a_sup) = [];

listmaxu = sign(diff(aux(1,:)));
ind_des_max = find(diff(listmaxu)<-1);
ind_des_max = ind_des_max+1;
maxdetectu = aux(2,ind_des_max);

ru = [1;maxdetectu(:);m+1];

end



%% Gaussian Scale Space
function [ fData ] = CreateGaussScaleSpace( data, width )
% space = CreateGaussScaleSpace( data, deriv, scales )
% Computes the Gaussian scale space of a signal
% Input:
%   data        A signal
%   width       the scale
% Output:
%   space       A scale space representation of the input data

% Find the gaussian kernel to convolve
g = GaussianKernel1D(width);
g = g/norm(g);

% we have to pad the data to avoid the derivative blowing up at the
% boundaries
padData = [data(1)*ones(3*width,1); data; data(end)*ones(3*width,1)];
fData = conv(padData, g);
% Resize the filtered data
offset = floor((length(fData) - length(data)) / 2);
fData = fData(offset+1:offset+length(data));
fData = fData/(range(fData(:)));

end


%% Gaussian function
function [ kernel ] = GaussianKernel1D( width  )
% GaussianKernel1D( width  ) Creates a centered Gaussian
% kernel of the first order and the specified scale. The
% width parameter specifies how many sigma from the center the kernel
% should extend.

sigma = max(round(width/3),1);
support = width;
range = -support:support;
center = range(support+1);
derivs(1:length(range)) = -((range-center)/(sigma^2));

kernel = (1/(sigma*sqrt(2*pi))) * exp(-((range-center).^2)/(2*sigma^2));
kernel = kernel .* derivs;

end
