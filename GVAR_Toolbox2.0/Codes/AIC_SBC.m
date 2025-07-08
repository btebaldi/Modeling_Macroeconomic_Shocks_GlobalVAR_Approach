function [logl aic sbc] = AIC_SBC(dep,X,bhat)

%**************************************************************************
% PURPOSE: Compute Akaike, Schwartz Bayesian statistics and log-likelihood 
%--------------------------------------------------------------------------
% INPUT:
% - dep = Tx1 series of dependent variable
% - X = Txm series of regressors
% - bhat = OLS coefficient vector (X'*X)\(X'*dep) 
%--------------------------------------------------------------------------
% OUTPUT:
% - logl = log-likelihood
% - aic = Akaike info criterion stat
% - sbc = Schwartz Bayesian info criterion stat
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************

T = rows(X);
s = rank(X);
res = dep-X*bhat;
Sigma=(1/T)*(res'*res);   % not adjusted
% Sigma=(1/(T-s))*(dep-xmat*A)'*(dep-xmat*A);  %degrees of freedom adjusted

% DdPS formulation
aic=(-T*(cols(dep)/2))*(1+log(2*pi))-(T/2)*log(det(Sigma))-cols(dep)*s;
sbc=(-T*(cols(dep)/2))*(1+log(2*pi))-(T/2)*log(det(Sigma))-(cols(dep)*(s/2))*log(T);
logl=(-T*(cols(dep)/2))*(1+log(2*pi))-(T/2)*log(det(Sigma));