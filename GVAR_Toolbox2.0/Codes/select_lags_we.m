function [logl_aic_sbc psc_degfrsc_Fcrit_Fsc] = select_lags_we(psc,endog_idv,ls,exog_idv,ln,exog2_idv,ec)

%**************************************************************************
% PURPOSE: Computes selection criteria for lag orders of weak exogeneity 
% regressions (AIC and SBC), together with F-statistics of residual serial 
% correlation tests.
%--------------------------------------------------------------------------
% INPUT: 
% - psc = max # of lagged residuals included in the regression
% - endog_idv: matrix, contains data of endogenous variables
% - ls: maximum lag order for the endogenous variables 
% - exog_idv: matrix, contains data of weakly exogenous variables 
% - ln: maximum lag order of weakly exogenous variables
% - exog2_idv: matrix, contains data of weakly exogenous variables 
% - ec: error correction term 
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
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************
     

dep = exog_idv -lagm(exog_idv,1);
dep = dep(2:end,:);

ch1 = endog_idv -lagm(endog_idv,1);
ch1 = ch1(2:end,:);
ch1block = [];
for i = 1:ls
    ch1block = [ch1block lagm(ch1,i)]; %#ok
end

ch2 = exog2_idv -lagm(exog2_idv,1);
ch2 = ch2(2:end,:);
ch2block = [];
for i = 1:ln
    ch2block = [ch2block lagm(ch2,i)]; %#ok
end


% now trim data according to the maximum lag order between domestic and
% foreign variables
lagtrim = max(ls,ln);
ch1block = trimr(ch1block,lagtrim,0);
ch2block = trimr(ch2block,lagtrim,0);
dep = trimr(dep,lagtrim,0);


% add ECM
ec = ec';
ec = trimr(ec,lagtrim,0);

% make sure that length of ECM is equal to length of dependent
% note: in general, length(ECM) >= length(dep)
discr = rows(dep) - rows(ec);
if discr > 0
    ch1block = trimr(ch1block,discr,0);
    ch2block = trimr(ch2block,discr,0);
    dep = trimr(dep,discr,0);
end

% add intercept
one = ones(rows(dep),1);

% block of regressors (includes EC terms)
X =[one ec ch1block ch2block];  

 % compute AIC, SBC
bhat = (X'*X)\(X'*dep);
[logl aic sbc] = AIC_SBC(dep,X,bhat);
logl_aic_sbc  = [logl aic sbc];


% F Tests for residual serial correlation
[degfrsc Fcrit Fsc] = Ftest_rsc(dep,X,psc);
psc_degfrsc_Fcrit_Fsc = [psc degfrsc Fcrit Fsc];


    
