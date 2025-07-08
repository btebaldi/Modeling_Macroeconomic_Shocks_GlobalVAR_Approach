function [kpsup kpmnsq ny rny qlr mw apw rqlr rmw rapw ...
          maxobs] = structural_stability_tests(cnum,cnames,endoglist,...
          endog,exog,varxlag,estcase,maxlag,ecm,ccut)

%**************************************************************************
% PURPOSE: do tests of structural stability of parameters
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
for n=1:cnum

    k = length(endoglist.(cnames{n}));
     
    y = endog.(cnames{n});
    x = exog.(cnames{n});
    
    lp = varxlag.(cnames{n})(1);
    lq = varxlag.(cnames{n})(2);
    
    Dy = y -lagm(y,1);    
    Dx = x -lagm(x,1);
    Z2= x -lagm(x,1);

    if lq~=1
        i=1;
        while i<=lq-1
            Z2 = [Z2 lagm(Dx,i)]; %#ok
            i = i+1;
        end
    end

    if lp~=1
        i=1;
        while i<=lp-1
            Z2 = [Z2 lagm(Dy,i)]; %#ok
            i=i+1;
        end
    end

    one = ones(rows(Z2),1);

    % trimming
    Z2 = trimr(Z2,maxlag,0);
    one = trimr(one,maxlag,0);


    if estcase.(cnames{n}) == 4 % case 4: Unrestricted intercepts; restricted trends
        Z2 = [one Z2]; %#ok
    elseif estcase.(cnames{n}) == 3 % case 3: Unrestricted intercepts in levels, no trends in Coint Space
        Z2 = [one Z2]; %#ok
    elseif estcase.(cnames{n}) == 2 % case 2: Restricted intercepts in Coint Space, no trends
        % do nothing (as the intercepts are restricted in the Coint Space)
    end

    
    Z2 = Z2';
    
    DY = trimr(Dy,maxlag,0);
    DX = [ecm.(cnames{n});Z2]';  % block of regressors


    s=1;
    while s<=k
        
        [kpsup_t kpmnsq_t] = kraplob(DY(:,s),DX);
        
        [ny_t rny_t] =nyblom(DY(:,s),DX);
        
        [qlr_t mw_t apw_t rqlr_t rmw_t rapw_t maxobs_t] = schow(DY(:,s),DX,ccut);
        
        endlist = endoglist.(cnames{n});
        
        kpsup.(endlist{s}).(cnames{n}) = kpsup_t;
        kpmnsq.(endlist{s}).(cnames{n}) = kpmnsq_t;
        ny.(endlist{s}).(cnames{n}) = ny_t;
        rny.(endlist{s}).(cnames{n}) = rny_t;
        qlr.(endlist{s}).(cnames{n}) = qlr_t;
        mw.(endlist{s}).(cnames{n}) = mw_t;
        apw.(endlist{s}).(cnames{n}) = apw_t;
        rqlr.(endlist{s}).(cnames{n}) = rqlr_t;
        rmw.(endlist{s}).(cnames{n}) = rmw_t;
        rapw.(endlist{s}).(cnames{n}) = rapw_t;
        maxobs.(endlist{s}).(cnames{n}) = maxobs_t;
        
        s=s+1;
    end
end