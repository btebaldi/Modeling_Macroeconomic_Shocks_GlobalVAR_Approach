
function [W chi2_pval] = jarquebera(x)

%**************************************************************************
% PURPOSE: Test for normality  (Jarque-Bera, Bowman-Shenton, 1977)
%--------------------------------------------------------------------------
% INPUT: 
% - x: T x 1 vector
%--------------------------------------------------------------------------
% OUTPUT:
% - W: Jarque-Bera statistic
% - chi2_pval: probability value of the computed J-B statistic
%--------------------------------------------------------------------------
% NOTE: it uses subroutine chis_cdf (Anders Holtsberg, 18-11-93)
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************


T=length(x);

mu=mean(x);
mu3=sum((x-mu).^3)/T;
mu4=sum((x-mu).^4)/T;
sig2=(var(x,1));
b1=(mu3/(sig2^(3/2)))^2;
b2=mu4/(sig2)^2;
W=T*( b1/6 + ((b2-3)^2)/24);

chi2_pval=1-chis_cdf(W,2);
