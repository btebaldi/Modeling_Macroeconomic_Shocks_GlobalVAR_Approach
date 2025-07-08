function  q = quantile(x,p)
%**************************************************************************
% PURPOSE: compute empirical quantile (percentile).
%--------------------------------------------------------------------------
% INPUT:
% x = matrix or vector of data
% p = chosen percentile
%--------------------------------------------------------------------------
% OUTPUT:
% q = empirical quantile


% USAGE:   q = quantile(x,p,method)
% where:   x = matrix or vector 
%          p = percent
%     method = 1,2,3
%    1. Interpolation so that F(X_(k)) == (k-0.5)/n (default)
%    2. Interpolation so that F(X_(k)) == k/(n+1).
%    3. The least number q such that at least a part p of x 
%      is less than or equal to q.     



if min(size(x)) == 1
   x = x(:);
   q = zeros(size(p));
else
   q = zeros(length(p),size(x,2));
end

x = sort(x); 
p = p(:);
n = size(x,1);

% -------------------------------------------------------------------------|
% The following method is used to match the results from Gauss.
% By Meiling He 

q(:) = x(int16(p*n),:);
