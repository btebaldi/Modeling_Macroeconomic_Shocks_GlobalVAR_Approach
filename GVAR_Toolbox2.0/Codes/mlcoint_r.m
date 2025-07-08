function [alpha Psi epsilon Omega ecm std ...
    logl aic sbc r2 rbar2 aics sbcs hcwstd ...
    nwcstd psc_degfrsc_Fcrit_Fsc] = mlcoint_r(beta_r,maxlag,estcase_idv,endog_idv,lp,exog_idv,lq,psc)

%**************************************************************************
% PURPOSE: It performs VECM/VECMX estimation under overidentifying
% restrictions
%--------------------------------------------------------------------------
% INPUT:
% - beta_r: overidentifying restrictions
% - maxlag: lag order of the GVAR
% - estcase_idv: treatment of deterministic components in the estimation:
% if =4, it is case IV (restricted trend in cointegration space,
% unrestricted intercept in levels); if =3, it is case III (no trend,
% unrestricted intercept in levels):  if =2, it is case II (no trend,
% restricted intercepts in cointegration space). Cases are from
% MacKinnon, Haug and Michelis (1999)
% - endog_idv: matrix containing domestic variables
% - lp: lag order of domestic variables
% - exog_idv: matrix containing foreign variables
% - lq: lag order of foreign variables
% - psc = max # of lagged residuals to use in serial correlation test
%--------------------------------------------------------------------------
% OUTPUT:
% - alpha: loading matrix coefficients estimates
% - Psi: short-run coefficients estimates
% - epsilon: residuals
% - Omega: estimated variance of VECM/VECMX model
% - ecm: error correction terms
% - std: standard errors of all coefficients estimates
% - logl: (restricted) log-likelihood of the VECM model
% - aic: Akaike statistic for the VECM model
% - sbc: Schwartz Bayesian statistic for the VECM model
% - r2: R-squared coefficient
% - rbar2: corrected R-squared coefficient
% - aics: Akaike statistics for each single equation of the VECM
% - sbcs: Schwartz Bayesian statistics for each single equation of the VECM
% - hcwstd: White's Heteroskedasticity robust standard errors
% - nwcstd: Newey-West HAC Consistent standard errors
% - psc_degfrsc_Fcrit_Fsc: serial correlation test on residuals
%--------------------------------------------------------------------------
% NOTES: if # input arguments is 5, VECM estimation is performed;
% if # input arguments is 7, VECMX estimation is performed.
% if # input arguments is 8, do VECMX and serial correlation test
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


if nargin == 5     % VECM estimation ( no exogenous variables )
    
    y = endog_idv;
    
    Z0 = y -lagm(y,1);
    
    Z1 = lagm(y,1);
    
    Z2 = [];
    
    if lp~=1
        i=1;
        while i<=lp-1
            Z2 = [Z2 lagm(Z0,i)]; %#ok
            i=i+1;
        end
    end
    
    % deterministic part
    trend = (1:rows(Z1))';
    one = ones(rows(Z1),1);
    
    % trimming
    Z0 = trimr(Z0,maxlag,0);
    Z1 = trimr(Z1,maxlag,0);
    if lp~=1
        Z2 = trimr(Z2,maxlag,0);
    end
    
    trend = trimr(trend,maxlag,0);
    one = trimr(one,maxlag,0);
    
elseif nargin == 7 || 8   % VECMX estimation ( exogenous variables included )
    
    
    y = endog_idv;
    x = exog_idv;
    
    Dy = y -lagm(y,1);
    Dx = x -lagm(x,1);
    
    Z0 = y -lagm(y,1);
    
    Z1=[lagm(y,1) lagm(x,1)];
    
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
    
    % deterministic part
    trend = (1:rows(Z1))';
    one = ones(rows(Z1),1);
    
    
    % trimming
    Z0 = trimr(Z0,maxlag,0);
    Z1 = trimr(Z1,maxlag,0);
    Z2 = trimr(Z2,maxlag,0);
    trend = trimr(trend,maxlag,0);
    one = trimr(one,maxlag,0);
    
end

if estcase_idv == 4 % case 4: Unrestricted intercepts; restricted trends
    Z1 = [trend-1 Z1];
    Z2 = [one Z2];
    
elseif estcase_idv == 3 % case 3: Unrestricted intercept in levels, no trend in Coint Space
    Z2 = [one Z2];
    
elseif estcase_idv == 2 % case 2: Restricted intercept in cointegration space, no trend
    Z1 = [one Z1];
    
end


% transpose series

Z0 = Z0';
Z1 = Z1';
Z2 = Z2';

[n T] = size(Z0);

% Compute product moment matrices Mij

if isempty(Z2)==0
    M02 = (1/T)*Z0*Z2';
    M12 = (1/T)*Z1*Z2';
    M22 = (1/T)*(Z2*Z2');
else
    M02 = [];
    M12 = [];
    M22 = [];
end

% Compute R(i,t) series
if isempty(Z2)==0;
    R0 = Z0-((M02/M22)*Z2);
    R1 = Z1-((M12/M22)*Z2);
else
    R0 = Z0;
    R1 = Z1;
end

% Compute product moment matrice Sij
S01 = (R0*R1')/T;
S11 = (R1*R1')/T;


% Impose betas:
beta = beta_r;

% Estimate parameters
alpha = (S01*beta)/(beta'*S11*beta);
Psi = (M02/M22)-((alpha*beta'*M12)/M22);

epsilon = Z0-(alpha*beta'*Z1)-(Psi*Z2);
Omega = (1/T)*(epsilon*epsilon');
ecm = beta'*Z1;

T=cols(epsilon); % # of observations

% Calculates standard errors
%***************************
residall = epsilon';

DX = [ecm;Z2]';  % block of regressors
s = cols(DX);

% non-robust standard errors
%**************************************

std = zeros(cols(residall),s);

for j=1:cols(residall)
    resid = residall(:,j);
    
    s2 = (1/(T-s))*(resid'*resid);
    varcov = s2*(eye(rows(DX'*DX))/(DX'*DX));
    
    std(j,:) = sqrt(diag(varcov));
end

% White's Heteroskedastic-robust Variance-Covariance Matrix
%*********************************************************

hcwstd = zeros(cols(residall),s);

for j=1:cols(residall)
    resid = residall(:,j);
    
    acc = 0;
    for i=1:rows(resid)
        acc = acc + (resid(i).*resid(i))*DX(i,:)'*DX(i,:);
    end
    
    HCWvarcov = (T/(T-s))*((DX'*DX)\acc)/(DX'*DX);
    
    hcwstd(j,:) = sqrt(diag(HCWvarcov));
    
end

% Newey-West HAC Consistent Variance-Covariance Matrix
%*********************************************************
nwcstd = zeros(cols(residall),s);

jw=1;
while jw<=cols(residall)
    seNW = neweywest(DX,residall(:,jw));
    nwcstd(jw,:)=seNW';
    jw=jw+1;
end


% VECMX models statistics
%**********************

% whole model statistics
% compute loglikelihood aic and sbc
logl=(-T*(n/2))*(1+log(2*pi))-(T/2)*log(det(Omega));
aic=(-T*(n/2))*(1+log(2*pi))-(T/2)*log(det(Omega))-n*s;
sbc=(-T*(n/2))*(1+log(2*pi))-(T/2)*log(det(Omega))-(n*(s/2))*log(T);


% single-equation statistics

% transpose dep variables
Z0 = Z0';

r2 = zeros(n,1);
rbar2 = zeros(n,1);
aics = zeros(n,1);
sbcs = zeros(n,1);

for i=1:n
    eps = epsilon(i,:)';
    ssr = eps'*eps;
    rsqr1 = ssr;
    ym = Z0(:,i) - ones(T,1)*mean(Z0(:,i));
    rsqr2 = ym'*ym;
    r2(i) = 1 - rsqr1/rsqr2;
    
    crsqr1 = rsqr1/(T-s);
    crsqr2 = rsqr2/(T-1);
    
    rbar2(i) = 1 - crsqr1/crsqr2;
    
    covmat = (1/T)*(eps'*eps);
    
    aics(i) =(-T/2)*(1+log(2*pi*diag(covmat)))-s;
    sbcs(i) =(-T/2)*(1+log(2*pi*diag(covmat)))-(s/2)*log(T);
end

if nargin == 8
    
    % F Tests for residual serial correlation
    %*************************
    dep = Z0;
    xmat = DX;
    Fsc = zeros(1,cols(dep));
    
    jsc = 1;
    while jsc <= cols(dep);        % making a regression for each dependent variable
        yidv=dep(:,jsc);
        xidv=xmat;
        
        Asc=(xidv'*xidv)\(xidv'*yidv);
        ressc=yidv-xidv*Asc;
        mx=eye(rows(yidv))-xidv*((xidv'*xidv)\xidv');
        wsc=zeros(rows(yidv),psc);
        
        isc=1;
        while isc <= psc;
            wsc(:,isc) = lagm(ressc,isc);
            isc = isc + 1;
        end
        
        chisqsc=rows(yidv)*((ressc'*wsc)*((wsc'*mx*wsc)\(wsc'*ressc))/(ressc'*ressc)); %/**compare with chi-square with p degrees of freedom**/
        Fsc(jsc)=((rows(yidv)-rank(xidv)-psc)/psc)*(chisqsc/(rows(yidv)-chisqsc)); %/** F(p,T-#regressors-p) **/
        degfrsc=(rows(yidv)-rank(xidv)-psc);
        jsc=jsc+1;
    end
    
    psc_degfrsc = [ psc degfrsc];
    
    % compute F-stat critical values at 95%
    if psc > 0
        Fcrit = fdis_inv(0.95,psc,degfrsc);
    else
        Fcrit = NaN;
    end
    
    psc_degfrsc_Fcrit_Fsc = [psc_degfrsc Fcrit Fsc  ];
    
else
    psc_degfrsc_Fcrit_Fsc = [];
    
end
