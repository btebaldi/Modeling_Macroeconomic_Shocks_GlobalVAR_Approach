function [a0_du a1_du Theta_du] =  vec2var_du(maxlag,gxvnum,varxlag_du,alpha_du,beta_du,Psi_du,estcase_du)
 

%**************************************************************************
% PURPOSE: Retrieve the VAR coefficients from the VEC dominant unit model
% estimates
%--------------------------------------------------------------------------
% INPUT:
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
% - gxvnum: number of global exogenous variables
% - varxlag_du: lag order of the dominant unit model
% - alpha_du: loadings, coefficients estimates
% - beta_du: cointegrating vector, coefficients estimates
% - Psi_du: short-run coefficients estimates
% - estcase_du: struct, contains treatment of deterministic components in the
% estimation: if =4, it is case IV (restricted trends in cointegration space,
% unrestricted intercepts in levels); if =3, it is case III (no trends in
% cointegration space, unrestricted intercepts in levels); if =2, it is case 
% II (no trends, restricted intercepts in cointegration space). Cases are from
% MacKinnon, Haug and Michelis (1999)
%--------------------------------------------------------------------------
% OUTPUT:
% - a0_du: intercept coefficients 
% - a1_du: linear trend slopes
% - Theta_du: coefficients of lagged endogenous variables 
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     

nall = gxvnum;
ndom = nall; % # of endogenous variables


lp = varxlag_du(1); % p lag of endogenous variables

%%%% extract intercept/trend coefficients
if estcase_du == 4
    a0_du = Psi_du(:,1); %  intercept (unrestricted estimates)
    Dcoeffs= Psi_du(:,2:end);  % leave intercept from Psi_du
    a1_du = alpha_du*beta_du(1,:)'; %  trend (restricted estimates)
    beta_du = beta_du(2:end,:); % exclude trend from beta_du (Cointegration vector)
elseif estcase_du == 3
    a0_du = Psi_du(:,1); %  intercept (unrestricted estimates)
    Dcoeffs= Psi_du(:,2:end);  % leave intercept from Psi_du
    a1_du = zeros(ndom,1); %  as no trend in case 3
elseif estcase_du == 2
    a0_du = alpha_du*beta_du(1,:)'; %  intercept (restricted estimates)
    Dcoeffs = Psi_du;
    a1_du = zeros(ndom,1); %  as no trend in case 2
    beta_du = beta_du(2:end,:); % exclude intercept from beta_du (Cointegration vector)
end
Dendcoeffs = Dcoeffs(:,end - (lp-1)*ndom +1:end);

%%%%
Pi = alpha_du*beta_du'; % Cointegration matrix


sumTheta_dusminuseye = Pi(:,1:ndom);
k=1;
Theta_ducoeffs_tmp = zeros(ndom,ndom,lp-1);
for j=1:lp-1
    Theta_ducoeffs_tmp(:,:,j) = Dendcoeffs(:,k:k+ndom-1);
    k=k+ndom;
end
if lp==1
    Theta_du(:,:,lp) = sumTheta_dusminuseye + eye(ndom);
elseif lp==2
    Theta_du(:,:,lp) = -Theta_ducoeffs_tmp(:,:,lp-1);
    Theta_du(:,:,1) = sumTheta_dusminuseye +Theta_ducoeffs_tmp(:,:,1) + eye(ndom);
else
    Theta_du(:,:,lp) = -Theta_ducoeffs_tmp(:,:,lp-1);
    for j=1:lp-2
        jj=lp-j;
        Theta_du(:,:,jj) = -(Theta_ducoeffs_tmp(:,:,jj-1)-(Theta_ducoeffs_tmp(:,:,jj)));
        
        
    end
    Theta_du(:,:,1) = sumTheta_dusminuseye + Theta_ducoeffs_tmp(:,:,1) + eye(ndom);
end

if lp < maxlag
    for j=1:maxlag-lp
        Theta_du(:,:,lp+j) = zeros(ndom,ndom);
    end
end






