function [delta_0 delta_1 H0 C zeta Sigma_zeta eta Sigma_eta eigens mods mlag] = solve_GVAR(maxlag,amaxlag,...
                          cnum,cnames,W,a0,a1,Theta,Lambda0,Lambda,Gamma0,Gamma,x,gxv,Wtilde,a0_du,a1_du,...
                          Theta_du,Lambda_du)

%**************************************************************************
% PURPOSE: It builds the GVAR model by stacking VARX country models and 
%          possibly the dominant unit model, then it calculates eigenvalues 
%--------------------------------------------------------------------------
% INPUT:
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
%           across country-specific models, it is the lag order of the GVAR 
%           if there is no dominant unit model
% - amaxlag: maximum lag order of endogenous and weakly exogenous variables
%            across country-specific models and dominant unit model, it is 
%            the lag order of the GVAR if there is a dominant unit model
% - cnum: number of country models
% - cnames: cell, list of countries' names (short names)
% - W: struct, contains country link matrices W_i 
% - a0: struct, contains intercept coefficients for each country model
% - a1: struct, contains linear trend slopes for each country model
% - Theta: struct, contains VARX coefficients of lagged endogenous 
%          variables for each country model
% - Lambda0: struct, contains VARX coefficients of contemporaneous weakly 
%            exogenous variables for each country model
% - Lambda: struct, contains VARX coefficients of lagged weakly exogenous
%           variables for each country model
% - Gamma0: struct, contains VARX coefficients of contemporaneous global
%           exogenous variables for each country model
% - Gamma: struct, contains VARX coefficients of lagged global exogenous
%          variables for each country model
% - x: matrix containing all domestic variables (global vector) 
% - gxv: T x gxvnum matrix containing all global exogenous variables 
%        (where gxvnum is the number of global exogenous variables)
% - Wtilde: fdvnum x K matrix containing weights used to construct the 
%           feedback variables of the dominant unit model, where fdvnum is
%           the maximum number of feedback variables across the augmented
%           regressions, while K is the number of domestic variables in the
%           GVAR
% - a0_du: gxvnum x 1 array which contains the intercept coefficients
%          of the dominant unit model
% - a1_du: gxvnum x 1 array which contains the trend coefficients
%          of the dominant unit model
% - Theta_du: struct, contains VARX coefficients of lagged global exogenous 
%             variables of the dominant unit model
% - Lambda_du: struct, contains VARX coefficients of lagged weakly exogenous 
%             variables (the feedback variables) of the dominant unit model 
%--------------------------------------------------------------------------
% OUTPUT:
% - delta_0: K x 1 array which contains stacked arrays of structure a0 
%            (it's the vector of intercepts of the GVAR)
% - delta_1: K x 1 array which contains stacked arrays of structure a1 
%            (it's the vector of trend slopes of the GVAR)
% - H0: K x K matrix containing impact coefficients of the GVAR
% - C: three-dim matrix (K x K x maxlag) containing estimated 
%      coefficients of the GVAR (i.e. of reduced form representation)
% - zeta: K x nobs matrix containing stacked country-specific residuals
% - Sigma_zeta: K x K varcov matrix obtained from residuals u
% - eta: K x nobs matrix containing residuals of the reduced-form GVAR 
% - Sigma_eta: K x K varcov matrix obtained from residuals eta
% - eigens: (K x maxlag) x 1 array of eigenvalues of the GVAR
% - mods: (K x maxlag) x 1 array of moduli of eigenvalues of the GVAR
% - mlag: lag order of the GVAR
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     

gxvnum = cols(gxv);

 
if gxvnum > 0 % there is a dominant unit model
    mlag = amaxlag;
else
    mlag = maxlag;
end
    

% Stacking VARX countries
%***********************
a_0 = []; % global intercept
a_1 = []; % global trend
for n=1:cnum
    a_0 = [a_0; a0.(cnames{n})]; %#ok
    a_1 = [a_1; a1.(cnames{n})]; %#ok
end
K = rows(a_0); % total endogenous variables in the system


G0=[];
for n=1:cnum
    if not(isempty(Lambda0.(cnames{n})))
        A.(cnames{n}) = [eye(rows(Lambda0.(cnames{n}))) -Lambda0.(cnames{n})];
    else
        A.(cnames{n}) = eye(rows(W.(cnames{n})));
    end
    G0=[G0; A.(cnames{n})*W.(cnames{n})]; %#ok
end


G = zeros(K,K,mlag);
for j=1:maxlag
    Gtemp = [];
    for n=1:cnum
        if not(isempty(Lambda0.(cnames{n})))
            B.(cnames{n})(:,:,j) = [Theta.(cnames{n})(:,:,j) Lambda.(cnames{n})(:,:,j)];
        else
            B.(cnames{n})(:,:,j) = Theta.(cnames{n})(:,:,j);
        end
        Gtemp = [Gtemp; B.(cnames{n})(:,:,j)*W.(cnames{n})]; %#ok     
    end
    G(:,:,j) = Gtemp;
end

if gxvnum > 0
    gexogvars = gxv';
    
    ts_gvar = cols(x);
    ts_du = cols(gexogvars);
    if ts_gvar > ts_du 
        x = x(:,1+ts_gvar-ts_du:end);
%     elseif ts_gvar < ts_du
%         gexogvars = gexogvars(:,1+ts_du-ts_gvar:end);
    end
    
    y = [x; gxv'];
    Ky = K + gxvnum;
else
    y = x;
    Ky = K;
end
yt = y';


trend_t =[0:1:cols(y)]'; %#ok
trend = trimr(trend_t,mlag,1);
nobst = rows(trimr(yt,mlag,0));


if gxvnum > 0 % there is a dominant unit model    
    % stack coefficients related to global exogenous variables
    J0=[];
    for n=1:cnum
        J0=[J0; Gamma0.(cnames{n})]; %#ok
    end
    H0 = [G0 -J0; zeros(gxvnum,K) eye(gxvnum)]; 
    
    J = [];
    Phi = [];
    H = [];
    for j=1:mlag
        % generate coefficient matrix of dominant unit model
        Phi(:,:,j) = Theta_du(:,:,j); %#ok
        % combine matrices of GVAR and dominant unit model
        if j<= maxlag
            Jtemp = [];
            for n=1:cnum
                Jtemp = [Jtemp; Gamma.(cnames{n})(:,:,j)]; %#ok
            end
            J(:,:,j) = Jtemp; %#ok
            if not(isempty(Wtilde))
            H(:,:,j) = [G(:,:,j) J(:,:,j); Lambda_du(:,:,j)*Wtilde Phi(:,:,j)]; %#ok    
            else
            H(:,:,j) = [G(:,:,j) J(:,:,j); zeros(gxvnum,K) Phi(:,:,j)]; %#ok
            end
        else % the GVAR model has lag order < than the dominant unit model
            if not(isempty(Wtilde))
            H(:,:,j) = [G(:,:,j) zeros(K,gxvnum); Lambda_du(:,:,j)*Wtilde Phi(:,:,j)]; %#ok    
            else
            H(:,:,j) = [G(:,:,j) zeros(K,gxvnum); zeros(gxvnum,K) Phi(:,:,j)]; %#ok
            end;
        end
    end
    
    h0 = [a_0; a0_du];
    h1 = [a_1; a1_du];
else
    H0 = G0;
    H = G;
    h0 = a_0;
    h1 = a_1;
end

acc = zeros(Ky,nobst);
for j=1:mlag
    acc = acc + H(:,:,j)*(trimr(lagm(yt,j),mlag,0))';
end
zeta = H0*(trimr(yt,mlag,0))' - h0*ones(1,nobst) - h1*trend' - acc;

 
Sigma_zeta = (zeta*zeta')/cols(zeta); 

% Building reduced form GVAR model parameters
%********************************************
delta_0 = H0\h0;
delta_1 = H0\h1;

C = [];
for j=1:mlag
    C(:,:,j) = H0\H(:,:,j); %#ok
end

% calculate companion form
Cbarup = [];
for j=1:mlag
    Cbarup = [Cbarup C(:,:,j)]; %#ok
end
Cbardown = [eye((K+gxvnum)*(mlag-1)) zeros((K+gxvnum)*(mlag-1),K+gxvnum)];
Cbar = [Cbarup; Cbardown];


eta = H0\zeta;
Sigma_eta = (eta*eta')/cols(eta);


% Calculating eigenvalues of the GVAR model

eigens = sort((eig(Cbar)),'descend');
mods = sort(abs(eig(Cbar)),'descend');

