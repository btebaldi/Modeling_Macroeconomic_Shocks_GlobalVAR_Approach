
function [lambda trace maxeig] = cointegration_test(maxlag,estcase,endog,lp,exog,lq)

%**************************************************************************
% PURPOSE: Calculate VAR or VARX eigenvalues, trace and maximum eigenvalue
% statistics in order to test the null hypothesis of r cointegrating
% relations
%--------------------------------------------------------------------------
% INPUT:
% - maxlag :    lag order of the GVAR
% - estcase:    treatment of deterministic components in the estimation: 
% if =4, it is case IV (restricted trend in cointegration space,
% unrestricted intercept in levels); if =3, it is case III (no trend in
% cointegration space, unrestricted intercept in levels); if = 2, it is case II
% (no trend, restricted intercept in cointegration space). Cases are from
% MacKinnon, Haug and Michelis (1999)
% - endog  :    matrix containing domestic variables
% - lp     :    lag order of domestic variables
% - exog   :    matrix containing foreign variables
% - lq     :    lag order of foreign variables
%--------------------------------------------------------------------------
% OUPUT:
% - lambda :    eigenvalues
% - trace  :    trace statistic
% - maxeig :    maximum eigenvalue statistic
%--------------------------------------------------------------------------
% NOTES: if # input arguments is 4, cointegration test for a VAR model is
% performed, elseif # input arguments is 6, it is performed for a VARX model
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************



if nargin == 4       % cointegration test for VAR model
    
    y = endog;
    
    Z0 = y - lagm(y,1);  % dependent variable
    
    Z1 = lagm(y,1);      % Ist block of regressors (cointegration space)
    
    Z2 = [];     % 2nd block of regressors
    if lp~=1
        i=1;
        while i<=lp-1
            Z2 = [Z2 lagm(Z0,i)]; %#ok
            i=i+1;
        end
    end
    
    % trimming in order to match # of observations
    Z0 = trimr(Z0,maxlag,0);
    Z1 = trimr(Z1,maxlag,0);
    if lp~=1
        Z2 = trimr(Z2,maxlag,0);
    end
    
elseif nargin == 6   % cointegration test for VARX model
    
    y = endog;
    x = exog;
    
    
    
    
    Z0=y - lagm(y,1); % dependent variable
    
    Z1=[lagm(y,1) lagm(x,1)]; % Ist block of regressors (cointegration space)
    
    Z2= x - lagm(x,1); % 2nd block of regressors
    
    if lp~=1
        i=1;
        while i<=lp-1
            Z2 = [Z2 lagm(Z0,i)]; %#ok
            i=i+1;
        end
    end
    
    Dx = x - lagm(x,1);
    
    if lq~=1
        i=1;
        while i<=lq-1
            Z2 = [Z2 lagm(Dx,i)]; %#ok
            i = i+1;
        end
    end
    
    % trimming in order to match # of observations
    Z0 = trimr(Z0,maxlag,0);
    Z1 = trimr(Z1,maxlag,0);
    Z2 = trimr(Z2,maxlag,0);
    
    
end

if estcase == 4
    % IV Case: Unrestricted intercepts; restricted trends
    
    trend = (1:rows(Z1))';
    Z1 = [trend-1 Z1];  % adding linear trend in cointegration space
    
    one = ones(rows(Z1),1);
    Z2 = [one Z2]; % adding intercept (unrestricted)
    
elseif estcase == 3
    % III Case: Unrestricted intercept in levels; no trend
    
    one = ones(rows(Z1),1);
    Z2 = [one Z2]; % adding intercept (unrestricted)
    
elseif estcase == 2
    % II Case: Restricted intercept in cointegration space; no trend
    
    one = ones(rows(Z1),1);
    Z1 = [one Z1]; % adding intercept (unrestricted)
end



% transpose series

Z0 = Z0';
Z1 = Z1';
Z2 = Z2';

[n T] = size(Z0);
nq1 = size(Z1,1);



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
S00 = (R0*R0')/T;
S01 = (R0*R1')/T;
S11 = (R1*R1')/T;
S10 = S01';


A = (S11\S10)*(S00\S01);   % estimation

[junk lambda] = eig(A); %#ok  % calculates eigenvalues and eigenvectors


lambda = diag(lambda);  % vectorize eigenvalues

lambda = sort(lambda,'descend'); % sort eigenvalues in descending order


% Remove zero eigenvalues for the case when the number of variables in Z1
% is greater than the number of variables in Z0 and adjust the eigenvectors
% accordingly

if n<nq1
    lambda = lambda(1:n);
end

% Trace statistic
trace = zeros(rows(Z0),1);
r=0;
tmp = zeros(rows(Z0),1);
for i= r+1:rows(Z0)
    tmp(i) = tmp(i) + (log(1 - lambda(i)));
    r = r+1;
end

r=0;
for i=1:rows(Z0)
    trace(i) = -T*(sum(tmp(1+r:end)));
    r=r+1;
end

% Maximum eigenvalue statistic
maxeig = zeros(rows(Z0),1);
for i = 1: rows(Z0)
    maxeig(i) = -T*log(1 - lambda(i));
end

