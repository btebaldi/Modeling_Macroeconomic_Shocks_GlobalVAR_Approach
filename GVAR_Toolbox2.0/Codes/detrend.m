
function resid = detrend(y,d)

%**************************************************************************
% PURPOSE: Given a time series y, it calculates residuals of a regression
% on deterministic components (intercept or intercept and linear trend)
%--------------------------------------------------------------------------
% INPUT: 
% - y : a time series column vector 
% - d : deterministic components index
%       if d = 0, performs a regression with intercept as regressor
%       if d = 1, performs a regression with intercept and linear trend as
%       regressors
%--------------------------------------------------------------------------
% OUTPUT:
% - resid: a column vector of residuals obtained from regression
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************



nobs = rows(y);

one = ones(nobs,1);
x = one;

if d == 1
    trend = (1:nobs)';
    x = [x trend];
end

b = (x'*x)\(x'*y);
resid = y - x*b;