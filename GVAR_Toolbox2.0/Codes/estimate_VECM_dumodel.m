function [a0_du a1_du Theta_du beta_du alpha_du Psi_du] = estimate_VECM_dumodel(maxlag,maxlag_du,varxlag_du,rank_du,estcase_du,esttype_du,gxv)

%**************************************************************************
% PURPOSE: estimating the dominant unit model
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************

gxvnum = cols(gxv);


if gxvnum == 1 % estimate AR(p) model
    
    if esttype_du == 0 % estimation in levels
        
        ptildel = varxlag_du(1);
        
        y = gxv(1+ptildel:end);
        if estcase_du == 0 % add intercept
            X = ones(rows(gxv),1);
        elseif estcase_du == 1 % add intercept and trend
            X = [ones(rows(gxv),1) [1:rows(gxv)]']; %#ok
        end
        i=1;
        while i<=ptildel
            X = [X lagm(gxv,i)]; %#ok
            i = i+1;
        end
        X = trimr(X,ptildel,0);
        
        bhat = (X'*X)\(X'*y);
        
        lp = ptildel;
        
        a0_du = bhat(1);
        if estcase_du == 0
            a1_du = 0; % no deterministic trend
            detind = 1;
        elseif estcase_du == 1
            a1_du = bhat(2);
            detind = 2;
        end
        Theta_du = zeros(1,1,lp);
        for i=1:maxlag
            if i<= lp
                Theta_du(:,:,i) = bhat(detind+i);
            else % if the lag order of the dominant unit model is smaller than the one of the gvar
                Theta_du(:,:,i) = 0;
            end
        end
        
        beta_du = [];
        alpha_du = [];
        Psi_du = [];
        
    elseif esttype_du == 1 % estimation in first differences
        
        ptildel = varxlag_du(1);

        y = gxv(1+ptildel:end) - gxv(ptildel:end-1);
        
        % add intercept by default
        X = ones(rows(gxv),1);
 
        i=1;
        while i<=ptildel
            X = [X lagm(y,i)]; %#ok
            i = i+1;
        end
        
        X = trimr(X,ptildel,0);
        y = trimr(y,ptildel,0);
        
        bhat = (X'*X)\(X'*y);
        
        
        lp = ptildel + 1;
        
        a0_du = bhat(1);
        a1_du = 0; % no deterministic trend
        Theta_du = zeros(1,1,maxlag);
        for i=1:maxlag
            if i == 1
                Theta_du(:,:,i) = 1 + bhat(1+i);     
            elseif i>1 && i<lp
                Theta_du(:,:,i) = bhat(1+i)-bhat(1+i-1);
            elseif i==lp
                Theta_du(:,:,i) = -bhat(1+i-1);
            else % if the lag order of the dominant unit model is smaller than the one of the gvar
                Theta_du(:,:,i) = 0;
            end
        end
        
        beta_du = [];
        alpha_du = [];
        Psi_du = [];
        
    end
    
else
    % VECM estimation
    [beta_du alpha_du Psi_du] = mlcoint(maxlag_du,rank_du,estcase_du,gxv,varxlag_du(1));
    
    % Retrieving the parameters of the VAR model from the VECM estimates
    %*******************************************************************
    ndom = rows(Psi_du); % # of endogenous variables
    lp = varxlag_du(1); % p lag of endogenous variables
  
    %- extract intercept coefficients
    a0_du = Psi_du(:,1); %  intercept (unrestricted estimates)
    Dcoeffs= Psi_du(:,2:end);  % leave intercept from Psi_du
    Dendcoeffs = Dcoeffs(:,end - (lp-1)*ndom +1:end);
    
    % - extract trend coefficients       
    if estcase_du == 4
    Alpha = alpha_du; % Loading matrix
    a1_du = Alpha*beta_du(1,:)'; %  trend (restricted estimates)
    Beta = beta_du(2:end,:); % exclude trend from beta (Cointegration vector)
    elseif estcase_du == 3
    Alpha = alpha_du; % Loading matrix
    a1_du = zeros(ndom,1); %  as no trend in case 3
    Beta = beta_du; % no need to exclude trend         
    end    
    
    % - recover coefficients related to lags of global exogeous variables
    Pi = Alpha*Beta'; % Cointegration matrix

    sumthetasminuseye = Pi(:,1:ndom);    
    k=1;
    thetacoeffs_tmp = zeros(ndom,ndom,lp-1);
    for j=1:lp-1
        thetacoeffs_tmp(:,:,j) = Dendcoeffs(:,k:k+ndom-1);
        k=k+ndom;
    end
    if lp==1
        Theta_du(:,:,lp) = sumthetasminuseye + eye(ndom);
    elseif lp==2
        Theta_du(:,:,lp) = -thetacoeffs_tmp(:,:,lp-1);
        Theta_du(:,:,1) = sumthetasminuseye +thetacoeffs_tmp(:,:,1) + eye(ndom);
    else
        Theta_du(:,:,lp) = -thetacoeffs_tmp(:,:,lp-1);    
    for j=1:lp-2
        jj=lp-j;
        Theta_du(:,:,jj) = -(thetacoeffs_tmp(:,:,jj-1)-(thetacoeffs_tmp(:,:,jj)));
    end
     Theta_du(:,:,1) = sumthetasminuseye + thetacoeffs_tmp(:,:,1) + eye(ndom);
    end
      
    if lp < maxlag
        for j=1:maxlag-lp
            Theta_du(:,:,lp+j) = zeros(ndom,ndom);
        end
    end      
end