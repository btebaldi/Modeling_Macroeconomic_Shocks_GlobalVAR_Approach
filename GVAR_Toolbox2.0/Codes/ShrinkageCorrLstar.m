function [Sigma_shrinkage lambda_star]=ShrinkageCorrLstar(Sig,T,lambda_param)

%**************************************************************************
% PURPOSE: Computes the lambda parameter for shrinkage of the
% variance-covariance matrix. 
%--------------------------------------------------------------------------
% INPUT:
% - Sig: unshrinked varcov matrix
% - T: number of time observations
% - lambda_param: shrinkage parameter, if empty, the routine computes it
% internally
%--------------------------------------------------------------------------
% OUPUT:
% - Sigma_shrinkage: shrinked variance-covariance matrix
% - lambda_star: shrinkage parameter
%--------------------------------------------------------------------------
% Vanessa Smith, 2014
% University of Cambridge, CFAP & CIMF, UK.
% lvs21@cam.ac.uk
%**************************************************************************

if not(isempty(lambda_param))
    lambda_star = lambda_param;
end

DmatR=diag((1./sqrt(diag(Sig))));
Rmat=DmatR*Sig*DmatR;

RmatM=Rmat+diag(NaN*ones(rows(Rmat),1));
vecRmatM=reshape(RmatM,rows(Rmat)*cols(Rmat),1);
vecRmatM(isnan(vecRmatM)) = [];
vecRmatSqM=vecRmatM.^2;

vecM=(vecRmatM.*(1-vecRmatSqM))/(2*T);
num=sum(vecRmatM.*(vecRmatM-vecM));
denom1=sum(((1-vecRmatSqM).^2)/T);
denom2=sum((vecRmatM-vecM).^2);
if isempty(lambda_param)
    lambda_star=1-(num/(denom1+denom2));
    if lambda_star<0
        lambda_star=0;
    elseif lambda_star>1
        lambda_star=1;
    end
end
Rmat_shrinkage=lambda_star*eye(rows(Sig))+(1-lambda_star)*Rmat;
DmatSigm=diag(sqrt(diag(Sig)));
Sigma_shrinkage=DmatSigm*Rmat_shrinkage*DmatSigm;
