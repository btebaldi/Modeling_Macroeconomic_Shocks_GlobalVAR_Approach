function [a0 a1 Theta Lambda0 Lambda] =  vecx2varx(maxlag,cnum,cnames,z,endoglist,varxlag,alpha,beta,Psi,estcase)

%**************************************************************************
% PURPOSE: Retrieve the parameters of each VARX country model from the VECMX
% estimates.
%--------------------------------------------------------------------------
% INPUT:
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
% - cnum: number of country models
% - cnames: cell, list of countries' names (short names)
% - z: struct, contains country series z_it = [domestic_it; foreign_it]
% - endoglist: struct, it contains lists of endogenous variables for each
% country
% - varxlag: struct, it contains lag orders for both endogenous and weakly
% exogenous variables for each country model
% - alpha: struct, contains loading matrix coefficients estimates
% - beta: struct, contains cointegrating vector coefficients estimates
% - Psi: struct, contains short-run coefficients estimates
% - estcase: struct, contains treatment of deterministic components in the
% estimation: if =4, it is case IV (restricted trends in cointegration space,
% unrestricted intercepts in levels); if =3, it is case III (no trends in
% cointegration space, unrestricted intercepts in levels); if =2, it is case 
% II (no trends, restricted intercepts in cointegration space). Cases are from
% MacKinnon, Haug and Michelis (1999)
%--------------------------------------------------------------------------
% OUTPUT:
% - a0: struct, contains intercept coefficients for each country model
% - a1: struct, contains linear trend slopes for each country model
% - Theta: struct, contains VARX coefficients of lagged endogenous 
% variables for each country model
% - Lambda0: struct, contains VARX coefficients of contemporaneous weakly 
% exogenous variables for each country model
% - Lambda: struct, contains VARX coefficients of lagged weakly exogenous
% variables for each country model
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     

for n=1:cnum
    
    nall = rows(z.(cnames{n})); % # of total variables for each country
    ndom = length(endoglist.(cnames{n})); % # of endogenous variables for each country
    nfor = nall - ndom; % # of weakly exogenous variables for each country

    lp = varxlag.(cnames{n})(1); % p lag of endogenous variables
    lq = varxlag.(cnames{n})(2); % q lag of exogenous variables


    Alpha = alpha.(cnames{n}); % Loading matrix
    %%%% extract intercept/trend coefficients       
    if estcase.(cnames{n}) == 4
    a0.(cnames{n}) = Psi.(cnames{n})(:,1); %  intercept (unrestricted estimates)
    Dcoeffs= Psi.(cnames{n})(:,2:end);  % leave intercept from Psi
    a1.(cnames{n}) = alpha.(cnames{n})*beta.(cnames{n})(1,:)'; %  trend (restricted estimates)
    Beta = beta.(cnames{n})(2:end,:); % exclude trend from beta (Cointegration vector)
    elseif estcase.(cnames{n}) == 3
    a0.(cnames{n}) = Psi.(cnames{n})(:,1); %  intercept (unrestricted estimates)
    Dcoeffs= Psi.(cnames{n})(:,2:end);  % leave intercept from Psi
    a1.(cnames{n}) = zeros(ndom,1); %  as no trend in case 3   
    Beta = beta.(cnames{n});
    elseif estcase.(cnames{n}) == 2
    a0.(cnames{n}) = alpha.(cnames{n})*beta.(cnames{n})(1,:)'; %  intercept (restricted estimates)
    Dcoeffs= Psi.(cnames{n});   
    a1.(cnames{n}) = zeros(ndom,1); %  as no trend in case 2
    Beta = beta.(cnames{n})(2:end,:); % exclude intercept from beta (Cointegration vector)        
    end
    Dendcoeffs = Dcoeffs(:,end - (lp-1)*ndom +1:end);
    Dforcoeffs = Dcoeffs(:,1:lq*nfor);
    
    
    %%%%
    Pi = Alpha*Beta'; % Cointegration matrix

    
    sumthetasminuseye = Pi(:,1:ndom);    
    k=1;
    thetacoeffs_tmp = zeros(ndom,ndom,lp-1);
    for j=1:lp-1
        thetacoeffs_tmp(:,:,j) = Dendcoeffs(:,k:k+ndom-1);
        k=k+ndom;
    end
    if lp==1
        Theta.(cnames{n})(:,:,lp) = sumthetasminuseye + eye(ndom);
    elseif lp==2
        Theta.(cnames{n})(:,:,lp) = -thetacoeffs_tmp(:,:,lp-1);
        Theta.(cnames{n})(:,:,1) = sumthetasminuseye +thetacoeffs_tmp(:,:,1) + eye(ndom);
    else
        Theta.(cnames{n})(:,:,lp) = -thetacoeffs_tmp(:,:,lp-1);    
    for j=1:lp-2
        jj=lp-j;
        Theta.(cnames{n})(:,:,jj) = -(thetacoeffs_tmp(:,:,jj-1)-(thetacoeffs_tmp(:,:,jj)));
    
    
    end
     Theta.(cnames{n})(:,:,1) = sumthetasminuseye + thetacoeffs_tmp(:,:,1) + eye(ndom);
    end
      
    if lp < maxlag
        for j=1:maxlag-lp
            Theta.(cnames{n})(:,:,lp+j) = zeros(ndom,ndom);
        end
    end  
    
    
    sumlambdas = Pi(:,ndom+1:end); % coefficient of x*_t-1 (foreign variables) in VECMX model    
    k=1;
    lambdacoeffs_tmp = zeros(ndom,nfor,lq);
    for j=1:lq
        lambdacoeffs_tmp(:,:,j) = Dforcoeffs(:,k:k+nfor-1);
        k=k+nfor;
    end
    Lambda0.(cnames{n}) = lambdacoeffs_tmp(:,:,1); 
    Lambda.(cnames{n})(:,:,lq) = -lambdacoeffs_tmp(:,:,lq);
    for j=1:lq-2
        jj=lq-j;
        Lambda.(cnames{n})(:,:,jj) = -lambdacoeffs_tmp(:,:,jj)+lambdacoeffs_tmp(:,:,jj+1);
    end
    
    acc = Lambda0.(cnames{n});
    for i=1:lq
        if not(i==1)
            acc = acc + Lambda.(cnames{n})(:,:,i);
        end
    end
    Lambda.(cnames{n})(:,:,1) = sumlambdas - acc;
    
    if lq < maxlag
        for j=1:maxlag-lq
            Lambda.(cnames{n})(:,:,lq+j) = zeros(ndom,nfor);
        end
    end
end





