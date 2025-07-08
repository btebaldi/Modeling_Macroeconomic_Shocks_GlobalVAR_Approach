function y = trimr(x,k1,k2)

%**************************************************************************
% PURPOSE: Trims observations from a time series array or matrix
%--------------------------------------------------------------------------
% INPUT:
% - x    : a time series array or matrix
% - k1   : # of initial observations to trim
% - k2   : # of final observations to trim
%--------------------------------------------------------------------------
% OUTPUT:
% - y    : array/matrix x after trimming
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
k = rows(x);

if (k1+k2) >= k
    error('Attempting to trim too much');
end

y = x(1+k1:k-k2,:);