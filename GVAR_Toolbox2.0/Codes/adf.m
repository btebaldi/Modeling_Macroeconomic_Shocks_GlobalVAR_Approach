
function out = adf(y,d,p)

%**************************************************************************
% PURPOSE: Performs Augmented Dickey-Fuller Test for Unit Root Presence
%--------------------------------------------------------------------------
% INPUT:
% - y: a nobs x 1 time series
% - d: indicates deterministic components to include in the regression
%      if d=0, it includes an intercept
%      if d=1, it includes an intercept and a linear trend
% - p: # of lagged changes to include in the regression
%--------------------------------------------------------------------------
% OUTPUT: out, containing:
% - dfstat: ADF statistic
% - aic:    Akaike Info Criterion value for the regression undertaken
% - sbc:    Schwartz Bayesian Info Criterion value for the regression
%           undertaken
% - d: indicates deterministic components to include in the regression
% - p: # of lagged changes to include in the regression
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

nobs = rows(y);    % # of observations

% Input control
if d<-1
    error('d cannot be set < -1');
elseif cols(y)>1
    error('ADF cannot handle a data matrix');
elseif (nobs - (2*p) + 1 <1)
    error('lag length is too large, negative degrees of freedom');
end

dep = diff(y,1); % dep: dependent variable
z = trimr(lagm(y,1),1,0); % z: block of regressors

if ( d > -1);
    one =ones(rows(z),1);
    if d==0   % Intercept only case
        z = [z one];
    elseif d==1   % Intercept and linear trend case
        trend =(1:rows(z))';
        z = [z one trend];
    end
end

% Creating lagged changes of y series to include in the regression
ch =  diff(y,1);
k=1;
while k <= p
    z = [z lagm(ch,k)]; %#ok
    k = k+1;
end

% Trimming
if p~=4
    z = trimr(z,k-1,0);
    dep = trimr(dep,k-1,0);
else
    z = trimr(z,k-2,0);
    dep = trimr(dep, k-2,0);
    z(1,end) = z(2,end);
end



% Estimation
b = (z'*z)\z'*dep; % coefficients
res = dep - z*b; % residuals

so  = (res'*res)/(rows(dep)-rank(z)); % Corrected estimated variance
soms = (res'*res)/(rows(dep)); % Not corrected estimated variance 

var_cov = so*(eye(rows(z'*z))/(z'*z));  % Variances-Covariances matrix

dfstat= b(1)/sqrt(var_cov(1,1)); % df statistic 

logl = -(rows(dep)/2)*((1+log(2*pi)) +log(soms));  % loglikelihood 

aic = -2*(logl/rows(dep)) + 2*(rank(z)/rows(dep)); % Akaike info crit
sbc = -2*(logl/rows(dep)) + rank(z)*log(rows(dep))/rows(dep); 
% Schwartz Bayesian info crit

% Storing output
out = [dfstat aic sbc logl]; 
