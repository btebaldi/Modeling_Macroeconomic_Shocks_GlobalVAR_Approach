function FEVDres = fevd(K,N,PHI,Sigma_u,G,eslct,sgirfflag,Sigma_zeta0)

%**************************************************************************
% PURPOSE: Computing forecast error variance decomposition 
%--------------------------------------------------------------------------
% INPUT:
% - K: total number of endogenous variables in the GVAR
% - N: forecast horizon
% - PHI: three-dim matrix (K x K x (N+1)) containing dynamic multipliers of
%   the GVAR
% - Sigma_u: K x K varcov matrix of the GVAR
% - G: K x K matrix containing impact coefficients of the GVAR
% - eslct: K x 1 array, selection vector for the shock to simulate
%--------------------------------------------------------------------------
% OUTPUT:
% - FEVDres: K x (N+1) matrix containing the generalized forecast error
%   variance decomposition
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************




if sgirfflag == 1 || sgirfflag == 2 % Structural GIRFs or OIRFs: orthogonalize
    
    P0_t = chol(Sigma_u(1:rows(Sigma_zeta0),1:rows(Sigma_zeta0)));
    P0 = P0_t';
    % so that we can orthogonalize residuals: v0 = P0*u0; ;
    
    %********
    P0_H0_tmp = eye(K);
    P0_H0_tmp(1:rows(P0),1:rows(P0)) = eye(rows(P0))/P0;
    P0_H0 = P0_H0_tmp;
    %********
    
    
    % apply Cholesky transformation
    Sigma_u = P0_H0*Sigma_u*P0_H0';
    % %
    % update G matrix
    G = P0_H0*G;
end


FEVDres = zeros(K,N+1);
vslct = eye(K);

eslct = abs(eslct);

invG = eye(rows(G))/G;
invGSigmau = G\Sigma_u;
%%%%%%

num = zeros(K,N+1);
den = zeros(K,N+1);

scale = 1./diag(Sigma_u);

n=1;
while n<=N+1
    for l=1:n
        acc1 = ((eslct'*PHI(:,:,l)*invGSigmau*vslct).^2)';
        num(:,n)= num(:,n) + acc1;
        acc2 = (eslct'*PHI(:,:,l)*invGSigmau*invG'*PHI(:,:,l)'*eslct);
        den(:,n)= den(:,n) + ones(K,1).*acc2;
    end
    FEVDres(:,n) = (scale.*num(:,n))./den(:,n);
    n=n+1;
end

