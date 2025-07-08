function  IRFres = irf(K,N,PHI,Sigma_u,G,eslct,sgirfflag,Sigma_zeta0)

%**************************************************************************
% PURPOSE: Computing impulse response functions 
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
% - IRFres: K x (N+1) matrix containing the generalized impulse response
%   functions
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

% preallocate matrix containing GIRFs
IRFres = zeros(K,N+1);

invGSigma_u = G\Sigma_u;

for i=1:N+1
    IRFres(:,i) = (PHI(:,:,i)*(invGSigma_u)*eslct)*(1/sqrt(eslct'*Sigma_u*eslct));
end

