function [degfrsc Fcrit Fsc] = Ftest_rsc(dep,X,psc)

%**************************************************************************
% PURPOSE: Perform F test for residual serial correlation
%--------------------------------------------------------------------------
% INPUT:
% - dep = Txk series of dependent variables
% - X = Txm series of regressors
% - psc = max # of lagged residuals to use in serial correlation test
%--------------------------------------------------------------------------
% OUTPUT:
% - degfrsc = number of degrees of freedom
% - Fcrit = critical value of the F test at the 95% significance level
% - Fsc = F test statistic
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************



T = rows(X);
Fsc = zeros(1,cols(dep));

jsc = 1;
while jsc <= cols(dep);        % making a regression for each dependent variable
    yidv=dep(:,jsc);

    Asc=(X'*X)\(X'*yidv);
    ressc=yidv-X*Asc;
    mx=eye(T)-X*((X'*X)\X');
    wsc=zeros(T,psc);

    isc=1;
    while isc <= psc;
        wsc(:,isc) = lagm(ressc,isc);
        isc = isc + 1;
    end

    chisqsc=T*((ressc'*wsc)*((wsc'*mx*wsc)\(wsc'*ressc))/(ressc'*ressc)); %/**compare with chi-square with p degrees of freedom**/
    Fsc(jsc)=((T-rank(X)-psc)/psc)*(chisqsc/(T-chisqsc)); %/** F(p,T-#regressors-p) **/
    degfrsc=(T-rank(X)-psc);
    jsc=jsc+1;
end

% compute F-stat critical values at 95%
if psc > 0
    
    try
        fdis_inv(0.95,psc,degfrsc);
    catch %#ok
        error('myApp:argChk', 'The F test critical values for testing serial correlation in the individual model residuals cannot be computed due to insufficient number of observations.')
    end
    
    Fcrit = fdis_inv(0.95,psc,degfrsc);
else
    Fcrit = NaN;
end