function B = lagm(A,p)

%**************************************************************************
% PURPOSE: given a time series matrix nobs x k, it yields the corresponding 
% nobs x k matrix of lagged series, where lag order is given by p
%--------------------------------------------------------------------------
% INPUT:
% - A: nobs x k matrix, where nobs is the number ot time observations
% - p: order of the lag operator
%--------------------------------------------------------------------------
% OUTPUT: 
% - B: a nobs x k matrix of lagged series 
%--------------------------------------------------------------------------
% EXAMPLE:
% given
% A= [1 2 3       
%     2 3 4
%     3 4 5
%     4 5 6]
% and p = 1
% it yields 
% B = [0 0 0
%      1 2 3
%      2 3 4
%      3 4 5]
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************


if nargin ==1
    p=1;
end

B = zeros(rows(A),cols(A));

for j=1:cols(A)
    for i=p+1:rows(A)
        B(i,j) = A(i-p,j);
    end
end

        