function [meanv medianv maxv minv sigma skew kurt]= dstats(x)

%**************************************************************************
% PURPOSE: computes descriptive statistics of a series
%--------------------------------------------------------------------------
% INPUT: 
% - x: T x 1 vector
%--------------------------------------------------------------------------
% OUTPUT:
% - meanv: scalar, Mean
% - medianv: scalar, Median
% - maxv: scalar, Maximum value
% - minv: scalar, Minimum value
% - sigma: scalar, standard deviation
% - skew: scalar, skewness coefficient
% - kurt: scalar, kurtosis coefficient
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

T=length(x);
meanv=mean(x);
medianv=median(x);
maxv=max(x);
minv=min(x);
sigma=std(x);
skew=sum(((x-meanv)/sigma).^3)/T;
kurt=sum(((x-meanv)/sigma).^4)/T;
