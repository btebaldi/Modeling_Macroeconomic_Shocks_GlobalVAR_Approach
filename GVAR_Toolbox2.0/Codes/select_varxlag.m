
function [logl_aic_sbc psc_degfrsc_Fcrit_Fsc] = select_varxlag(maxlag,psc,endog,lp,exog,lq)

%**************************************************************************
% PURPOSE: Creates VARX models, then calculates AIC and SBC values, and 
%          performs F-tests for residual serial correlation
%--------------------------------------------------------------------------
% INPUT: 
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
%           for each country
% - psc = max # of lagged residuals included in the regression
% - endog: structure, contains data of endogenous variables of each country
% - lp: lag order of endogenous variables in the VARX
% - exog: structure, contains data of weakly exogenous variables of each country
% - lq: lag order of weakly exogenous variables in the VARX
%--------------------------------------------------------------------------
% OUTPUT:
% - logl_aic_sbc : array, contains loglikelihood, AIC and SBC statistics
% - psc_degfrsc_Fcrit_Fsc: array, contains 
%      - psc = # of lagged residuals included in the regression
%      - degfrsc = # degrees of freedom (e.g # of observations - #
%        of regressors - # of lagged residuals in the regression)
%      - Fcrit = F-critical value at 95% 
%      - Fsc = F-stat of the test for serial correlation
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     


if (lp<0) || (lq<0)
    error(' negative lag order ');
end


nobs = rows(endog);

one = ones(nobs,1);
trend = (1:nobs)';
xmat = [one trend];

if lp~=0
    i=1;
    while i<=lp
        xmat = [ xmat lagm(endog,i)]; %#ok
        i=i+1;
    end
end

xmat = [xmat exog];

if (lq~=0)
    i=1;
    while i <= lq
        xmat = [xmat lagm(exog,i)]; %#ok
        i=i+1;
    end
end


% trimming
dep = trimr(endog,maxlag,0);
xmat = trimr(xmat,maxlag,0);

% OLS
A=(xmat'*xmat)\(xmat'*dep);

% AIC - BIC
[logl aic sbc] = AIC_SBC(dep,xmat,A);
logl_aic_sbc  = [logl aic sbc];


% F Tests for residual serial correlation
[degfrsc Fcrit Fsc] = Ftest_rsc(dep,xmat,psc);
psc_degfrsc_Fcrit_Fsc = [psc degfrsc Fcrit Fsc];

