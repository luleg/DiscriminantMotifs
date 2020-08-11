function [r,c,res,MVP] = bnewtns(A,p,q,tol,delta,fl,Delta)

% BNEWTNS A balancing algorithm for nonsymmetric matrices
%
% [R, C] = BNEWTNS(A) attempts to find vectors r and c such that
% diag(R)*A*diag(C) is close to bistochastic.
% A must be nonnegative.

% [R, C] = BNEWTNS(A, P, Q) attempts to find vectors r and c such that
% diag(R)*A*diag(C) has row sums close to P and column sumes close to Q.
%
% [R, C] = BNEWTNS(A,P,Q,TOL) specifies the tolerance of the method.
% If TOL = [] then BNEWT uses the default, 1e-6.
%
% [R, C] = BNEWTNS(A,P,Q,TOL,DEL) determines how close we are
% willing to allow our balancing vectors can get to the edge of the
% positive cone. We use a relative measure on the size of elements. The
% default value for DEL is 0.1
%
% [R, C] = BNEWTNS(A,P,Q,TOL,DEL,FL) will output intermediate convergence
% statistics if FL = 1. The default value for FL is 1.
%
%
% [R, C, RES] = BNEWTNS(A,P,Q,TOL,DEL,FL) returns the residual error, too.
[m,n]=size(A);
if nargin == 1,p  = ones(m,1); q = ones(n,1); end
e = ones(m+n,1); res=[]; h=[p;q];
if nargin < 7, Delta = 3; end
if nargin < 6,  fl = 0; end
if nargin < 5,  delta= 0.1; end
if nargin < 4,  tol = 1e-6; end


g=0.9; etamax = 0.1; % Parameters used in determining inner iteration
% stopping criterion. Default g=0.9; eta = 0.1;
eta = etamax; stop_tol = tol*.5;  %Relative tolerance?

x = e; rt = tol^2;
v = x.*[A*x(m+1:end); A'*x(1:m)]; rk = h - v;
rho_km1 = rk'*rk; rout = rho_km1; rold = rout;

MVP = 0;  % We'll count matrix vector products.
i = 0; % Outer iteration count.

if fl == 1, fprintf('it    in. it    res\n'), end
while rout > rt && MVP < 25e3 % Outer iteration
    i = i + 1; k = 0; y = e;
    innertol = max([eta^2*rout,rt]);
    while  rho_km1 > innertol %Inner iteration by CG
        k = k + 1;
        if k == 1
            Z = rk./v;  p=Z; rho_km1 = rk'*Z;
        else
            beta=rho_km1/rho_km2;
            p=Z + beta*p;
        end
        tx = x.*p; w = x.*[A*tx(m+1:end); A'*tx(1:m)] + v.*p;
        alpha = rho_km1/(p'*w);
        ap = alpha*p;
        % We want to avoid -ve components. To achieve this we limit how
        % much closer we can get to the boundary in each inner iteration.
        ynew = y + ap;
        if min(ynew) <= delta
            if delta == 0, break, end
            ind = find(ap < 0);
            gamma = min((delta  - y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        if max(ynew) >= Delta
            ind = find(ynew > Delta);
            gamma = min((Delta-y(ind))./ap(ind));
            y = y + gamma*ap;
            break
        end
        y = ynew;
        rk = rk - alpha*w;  rho_km2 = rho_km1;
        Z = rk./v;  rho_km1 = rk'*Z;
    end
    x = x.*y; v = x.*[A*x(m+1:end); A'*x(1:m)];
    rk = h - v; rho_km1 = rk'*rk; rout = rho_km1;
    
    MVP = MVP + k + 1;
    
    % Update inner iteration stopping criterion.
    rat = rout/rold;  rold = rout; res_norm = sqrt(rout);
    eta_o = eta;  eta = g*rat;
    if g*eta_o^2 > 0.1
        eta = max([eta,g*eta_o^2]);
    end
    eta = min([eta,etamax]);
    eta = max([eta,stop_tol/ res_norm]);
    
    if fl == 1
        fprintf('%3d %6d   %.3e %.3e %.3e %.3e\n', i,k, res_norm,min(y),max(y),min(x));
        res=[res; res_norm];
    end
end

% fprintf('Matrix-vector products = %6d\n', MVP)
r = x(1:m); c = x(m+1:end);

end