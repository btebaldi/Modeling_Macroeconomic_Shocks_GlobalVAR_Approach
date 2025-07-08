function wetest_stat = test_weakexogeneity(cnames,cnum,endog,exog,varxlag_we,exog_we,ecm,rank)

                  
%**************************************************************************
% PURPOSE: Performing weak exogeneity test 
%--------------------------------------------------------------------------
% INPUT: 
% - cnames: cell, list of countries' names (short names)
% - cnum: number of country models
% - endog: structure, contains data of endogenous variables of each country
% - exog: structure, contains data of weakly exogenous variables of each 
% country
% - varxlag_we: structure, contains lag orders for weak exogeneity test
% regressions for each country
% - exog_we: structure, contains data of weakly exogenous variables for 
% weak exogeneity test regressions for each country
% - ecm: structure, contains error correction terms of each country
% - rank: structure, contains ranks of each country
%--------------------------------------------------------------------------
% OUTPUT:
% - wetest_stat: structure: it contains F-stats of weak exogeneity test,
% corresponding degrees of freedom and the critical values at 95%  
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     


for n=1:cnum

    endog_idv = endog.(cnames{n});
    exog_idv = exog.(cnames{n});
    exog2_idv = exog_we.(cnames{n});
    
    ls = varxlag_we.(cnames{n})(1);
    ln = varxlag_we.(cnames{n})(2);
    
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
    ec = ecm.(cnames{n})';
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
    
    % F Tests
    Fstat = zeros(1,cols(dep));
    
    for j=1:cols(dep)       % making a regression for each dependent variable
        Y=dep(:,j);
             
        % unrestricted regression
        X=[one ec ch1block ch2block];  % block of regressors (includes EC terms)
        A=(X'*X)\(X'*Y);
        u=Y-X*A;
        
        % restricted regression ( ECM coefficients imposed = 0)
        Xr = [one ch1block ch2block];  % block of regressors (excludes EC terms)
        Ar = (Xr'*Xr)\(Xr'*Y);
        ur = Y-Xr*Ar;
        
        % compute F-statistic:
        restr = rank.(cnames{n}); %number of restrictions
        k = cols(X); %number of regressors
        T = rows(Y);  % number of observations
        
        Fstat(j)=((ur'*ur - u'*u)/restr)/((u'*u)/(T-k));
        degfr=T-k;
        
    end
    
    % compute F-stat critical values at 95%
    if restr > 0
        fcrit = fdis_inv(0.95,restr,degfr);
    else
        fcrit = NaN;
    end
    
    wetest_stat.(cnames{n}) = [restr degfr fcrit Fstat];
end

