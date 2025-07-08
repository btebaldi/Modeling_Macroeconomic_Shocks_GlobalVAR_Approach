function [ptildel_found qtildel_found tab a0_du a1_du Theta_du Lambda_du amaxlag X dep] = augmentedregression(isfeedback,lagselect_du2,info2,max_ptildel,max_qtildel,maxlag,...
         gxv,gxvnames,exog_du,exognames_du,feedbacks_flagmat,ecm_du,alpha_du,beta_du,estcase_du,esttype_du2,estcase_du2,psc_du2)


%**************************************************************************
% PURPOSE: Estimate the augmented regression (II stage) of the dominant 
%          unit model
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************



gxvnum = cols(gxv);
ncells = 100;
ptildel_found = [];
qtildel_found = [];



% for printing results, in case lag selection == 1
vlabel = [];
vres = [];
flabel = [];
fres = [];

if isfeedback == 0
    max_qtildel = ones(gxvnum,1); % define max_qtildel, even though irrelevant
    % as there are no feedbacks
end


Feqspec = [];
for g = 1:gxvnum
    if lagselect_du2 == 1 % do optimal selection of lag orders
        ini_ptildel = 1;
        ini_qtildel = 1;
    elseif lagselect_du2 == 0 % don't do lag order selection
        ini_ptildel = max_ptildel(g);
        if not(isempty(exog_du.(gxvnames{g})))
        ini_qtildel = max_qtildel(g);
        else
        ini_qtildel = 1;
        end
    end
    
    resmatrix = [];
    for j1 = ini_ptildel:max_ptildel(g)
        ptildel = j1;

        for j2 = ini_qtildel:max_qtildel(g)
            qtildel = j2;
            
            % model specification
            %**************************************************************
            
            % specification of the dependent variable
            if gxvnum == 1 % onevariate model
                if esttype_du2 == 0 % estimation in levels
                    dep = gxv;
                elseif esttype_du2 == 1 % % estimation in first differences
                    dep = gxv-lagm(gxv);
                    dep = dep(2:end); % trim first obs
                end
                % in this case the set of all global exogenous variables
                % corresponds to the dependent variable
                gxvset = dep;
            elseif gxvnum > 1 % multivariate model
                dep = gxv(:,g)-lagm(gxv(:,g)); 
                dep = dep(2:end); % trim first obs
                gxvset = gxv-lagm(gxv);
                gxvset = gxvset(2:end,:); % trim first obs
            end
            
            % specification of the regressors
            if gxvnum == 1 % onevariate model 
                if estcase_du2 == 0 % add intercept
                    X = ones(rows(dep),1);
                elseif estcase_du2 == 1 % add intercept and trend
                    X = [ones(rows(dep),1) [1:rows(dep)]']; %#ok
                end
                % check if feedbacks are included
                if not(isempty(exog_du.(gxvnames{1}))) % there are feedbacks
                    if esttype_du2 == 0 % estimation in levels
                        % feedbacks enter in levels
                        fback = exog_du.(gxvnames{1});
                    elseif esttype_du2 == 1 % % estimation in first differences
                        % feedbacks enter in first differences
                        fback = exog_du.(gxvnames{1})-lagm(exog_du.(gxvnames{1}),1);
                        fback = fback(2:end,:); % trim first obs
                    end
                end
            elseif gxvnum > 1 % multivariate model
                % add intercept by default
                X = ones(rows(dep),1);
                
                % check if feedbacks are included
                if not(isempty(exog_du.(gxvnames{g}))) % there are feedbacks
                    % feedbacks enter in first differences
                    fback = exog_du.(gxvnames{g})-lagm(exog_du.(gxvnames{g}),1);
                    fback = fback(2:end,:); % trim first obs
                end
            end
            
            % add the lags of all set of global exogenous variables as regressors (and if chosen, lags of feedbacks)
            i=1;
            while i<=ptildel
                X = [X lagm(gxvset,i)]; %#ok
                i = i+1;
            end
            
            if not(isempty(exog_du.(gxvnames{g})))
                i=1;
                while i<=qtildel
                    X = [X lagm(fback,i)]; %#ok
                    i = i+1;
                end
            end   
            
            if (gxvnum > 1) && (not(isempty(ecm_du))) % add the error correction term
                ecmterm = ecm_du';
                % now add the ecm term
                T_ecm = rows(ecmterm);
                T_other = rows(X);
                if T_other > T_ecm
                    dep = trimr(dep,T_other-T_ecm,0);
                    X = trimr(X,T_other-T_ecm,0);
                elseif T_other < T_ecm
                    ecmterm = trimr(ecmterm,T_ecm-T_other,0);
                end
                X = [X ecmterm]; %#ok
            end
            
            % trimming
            if not(isempty(exog_du.(gxvnames{g})))
                maxtrim = max(ptildel,qtildel);
            else
                maxtrim = ptildel;
            end
            X = trimr(X,maxtrim,0); 
            dep = trimr(dep,maxtrim,0);
          
            
            % do OLS
            bhat = (X'*X)\(X'*dep);
            
            % compute AIC, SBC
            [logl aic sbc] = AIC_SBC(dep,X,bhat);
            
            % F Test for residual serial correlation
            [degfrsc Fcrit Fsc] = Ftest_rsc(dep,X,psc_du2);
            
            resmatrix = [resmatrix; j1 j2 aic sbc psc_du2 degfrsc Fcrit Fsc]; %#ok
            
            % for printing results
            %*********************
            vlabel = [vlabel; gxvnames(g)]; %#ok  
            if not(isempty(exog_du.(gxvnames{g})))
                vres = [vres; j1 j2 aic sbc logl]; %#ok
            else
                vres = [vres; j1 NaN aic sbc logl]; %#ok
            end
            flabel_tmp = sprintf('F(%d,%d)', psc_du2,degfrsc);
            flabel_tmp = {flabel_tmp};            
            flabel = [flabel; flabel_tmp]; %#ok
            fres = [fres; Fcrit Fsc]; %#ok
        end
    end
    

    if lagselect_du2 == 1 % recover optimally chosen lag orders
        if info2 == 2 % aic
            [junk bestfit_idx] = max(resmatrix(:,3)); %#ok
        elseif info2 == 3 % sbc
            [junk bestfit_idx] = max(resmatrix(:,4)); %#ok
        end
        ptildel_found = [ptildel_found; resmatrix(bestfit_idx,1)]; %#ok
        qtildel_found = [qtildel_found; resmatrix(bestfit_idx,2)]; %#ok
    elseif lagselect_du2 == 0 % store the OLS estimates
        
        olsest.(gxvnames{g}) = bhat;
        
        % just for printing purposes define
        est_det.(gxvnames{g}) = bhat(1);
        if (gxvnum > 1) && (not(isempty(ecm_du)))
            est_ecm.(gxvnames{g}) = bhat(end+1-cols(ecmterm):end);
            est_oth.(gxvnames{g}) = bhat(2:end-cols(ecmterm));
        else
            est_ecm.(gxvnames{g}) = [];
            est_oth.(gxvnames{g}) = bhat(2:end);
        end

        % store the F-stats
        Fcommon = [psc_du2 degfrsc Fcrit];
        Feqspec = [Feqspec Fsc]; %#ok
    end
end

% determine the lag orders of the dominant unit model

% - for the global variables
if gxvnum == 1 && esttype_du2 == 0 % univariate model estimated in levels
    ptilde = max(max_ptildel);
else
    ptilde = max(max_ptildel) + 1;
end


% for the feedbacks, if included
if not(isempty(exog_du))
    if gxvnum == 1 && esttype_du2 == 0 % univariate model estimated in levels
        qtilde = max(max_qtildel);
    else
        qtilde = max(max_qtildel) + 1;
    end
else
    qtilde = 0;
end

% determine the lag order of the "augmented" gvar, as
% max(maxlag,ptilde,qtilde)
amaxlag = max([maxlag ptilde qtilde]);


if lagselect_du2 == 0 % recover the estimates of the dominant unit model
    tab = ['OLS Estimates of Augmented Regressions' num2cell(NaN(1,ncells-1))];
    tab = [tab; num2cell(NaN(1,ncells))];
    
    if gxvnum == 1 % onevariate model, we have just a set of OLS estimates, bhat
        
        ptildel = max_ptildel(1);
        qtildel = max_qtildel(1);
        
        if esttype_du2 == 0 % estimation in levels
            
            a0_du = bhat(1);
            if estcase_du2 == 0
                a1_du = 0; % no deterministic trend
                detind = 1;
            elseif estcase_du2 == 1
                a1_du = bhat(2);
                detind = 2;
            end                   
            
            % feed coefficients of lagged dependent variable
            Theta_du = zeros(1,1,amaxlag);
            for i=1:amaxlag
                if i<= ptilde
                    Theta_du(:,:,i) = bhat(detind+i);
                else % if the lag order of the dominant unit model is smaller than the one of the augmented gvar
                    Theta_du(:,:,i) = 0;
                end
            end
            
            
            % feed coefficients of feedbacks, if included
            xtildelab = []; 
            if isfeedback == 1
                maxfeednum = 0; feedbacks_flagmat_redux = [];
                for j=1:cols(feedbacks_flagmat)
                    pointer = sum(feedbacks_flagmat(:,j));
                    if pointer > 0
                        maxfeednum = maxfeednum + 1;
                        feedbacks_flagmat_redux = [feedbacks_flagmat_redux feedbacks_flagmat(:,j)]; %#ok
                    end
                end
                Lambda_du = zeros(1,maxfeednum,amaxlag);
                
                k=1;
                for i=1:amaxlag
                    if i<= qtilde
                        aj=1;
                        for j=1:cols(feedbacks_flagmat_redux)
                            if feedbacks_flagmat_redux(g,j) == 1
                                Lambda_du(1,j,i) = bhat(detind+ptilde+k);
                                k=k+1;

                                xtildelab_tmp = sprintf('%s論%d',exognames_du.(gxvnames{1}){aj},i);
                                xtildelab_tmp = {xtildelab_tmp};
                                xtildelab = [xtildelab xtildelab_tmp]; %#ok
                                
                                aj=aj+1; 
                            end
                        end   
                    end
                end
            else
                Lambda_du = [];
            end
              
        elseif esttype_du2 == 1 % estimation in first differences
            
            a0_du = bhat(1);
            a1_du = 0; % no deterministic trend
            detind = 1;
            % feed coefficients of lagged dependent variable
            for i=1:amaxlag
                if i == 1
                    Theta_du(:,:,i) = 1 + bhat(detind+i); %#ok
                elseif i>1 && i<ptilde
                    Theta_du(:,:,i) = bhat(detind+i)-bhat(detind+i-1); %#ok
                elseif i==ptilde
                    Theta_du(:,:,i) = -bhat(detind+i-1); %#ok
                else % if the lag order of the dominant unit model is smaller than the one of the gvar
                    Theta_du(:,:,i) = 0; %#ok
                end
            end
            
            
            % feed coefficients of feedbacks, if included
            if isfeedback == 1
                maxfeednum = 0; feedbacks_flagmat_redux = [];
                for j=1:cols(feedbacks_flagmat)
                    pointer = sum(feedbacks_flagmat(:,j));
                    if pointer > 0
                        maxfeednum = maxfeednum + 1;
                        feedbacks_flagmat_redux = [feedbacks_flagmat_redux feedbacks_flagmat(:,j)]; %#ok
                    end
                end
                Lambda_du = zeros(1,maxfeednum,amaxlag);
                
                k=1;
                for i=1:amaxlag
                    if i == 1
                        for j=1:cols(feedbacks_flagmat_redux)
                            if feedbacks_flagmat_redux(g,j) == 1
                                Lambda_du(1,j,i) = bhat(detind+ptilde-1+k);
                                k=k+1;
                            end
                        end
                    elseif i>1 && i<qtilde
                        for j=1:cols(feedbacks_flagmat_redux)
                            if feedbacks_flagmat_redux(g,j) == 1
                                Lambda_du(1,j,i) = bhat(detind+ptilde-1+k)-bhat(detind+ptilde-1+k-maxfeednum);
                                k=k+1;
                            end
                        end                        
                    elseif i== qtilde
                        k = k-maxfeednum;
                        for j=1:cols(feedbacks_flagmat_redux)
                            if feedbacks_flagmat_redux(g,j) == 1
                                Lambda_du(1,j,i) = -bhat(detind+ptilde-1+k);
                                k=k+1;
                            end
                        end                            
                    end
                end
            else
                Lambda_du = [];
            end                
        end
        
        % prepare output matrices
        %---------------------------
        
        % labels of dependent variable and deterministic components
        if esttype_du2 == 0
        deplab_tmp = sprintf('%s',gxvnames{1});    
        elseif esttype_du2 == 1
        deplab_tmp = sprintf('d%s',gxvnames{1});
        end
        deplab_tmp = {deplab_tmp};
        deplab = deplab_tmp;
        
        detlab = {'Intercept'};
        if esttype_du2 == 0 && estcase_du2 == 1
            detlab = [detlab 'Trend'];
        end
        
        
        % feed coefficients of lagged dependent variable
        xlab = [];
        for i=1:ptildel
            if esttype_du2 == 0
            xlab_tmp = sprintf('%s_%d',gxvnames{1},i);
            elseif esttype_du2 == 1
            xlab_tmp = sprintf('d%s_%d',gxvnames{g},i);
            end
            xlab_tmp = {xlab_tmp};
            xlab = [xlab xlab_tmp]; %#ok
        end
        
        xtildelab = [];
        if not(isempty(exog_du.(gxvnames{1})))
            for q = 1:qtildel
                for i=1:cols(exog_du.(gxvnames{1}))
                    if esttype_du2 == 0
                    xtildelab_tmp = sprintf('%s論%d',exognames_du.(gxvnames{1}){i},q);    
                    elseif esttype_du2 == 1
                    xtildelab_tmp = sprintf('d%s論%d',exognames_du.(gxvnames{1}){i},q);
                    end
                    xtildelab_tmp = {xtildelab_tmp};
                    xtildelab = [xtildelab xtildelab_tmp]; %#ok
                end
            end
        end
        
        
        hlab = [deplab detlab xlab xtildelab]; 
        hres = [NaN bhat'];
        
        tab = [tab; hlab num2cell(NaN(1,ncells-length(hlab))); num2cell(hres) num2cell(NaN(1,ncells-length(hres)))];
        
    elseif gxvnum > 1  % multivariate model, we have gxvnum sets of OLS estimates
        
        a0_du = zeros(gxvnum,1); % initialize matrix of intercepts
        
        if (not(isempty(ecm_du)))
            %%%% extract intercept/trend coefficients
            if estcase_du == 4
                a1_du = alpha_du*beta_du(1,:)'; %  trend (restricted estimates)
                Beta = beta_du(2:end,:); % exclude trend from beta (Cointegration vector)
            elseif estcase_du == 3
                a1_du = zeros(gxvnum,1); %  as no trend in case 3
                Beta = beta_du; % no need to exclude trend
            elseif estcase_du == 2
                a0_du = alpha_du*beta_du(1,:)'; %  intercept (restricted estimates)
                a1_du = zeros(gxvnum,1); %  as no trend in case 2
                Beta = beta_du(2:end,:); % exclude intercept
            end
        else
            a1_du = zeros(gxvnum,1); %  as no trend in case 2
        end
        
        % OLS coefficients related to the error correction terms
        rank_du = cols(alpha_du);
        %delta_du = zeros(gxvnum,rank_du); 
        
        % OLS coefficients related to the lagged changes of 
        % global exogenous variables 
        Gamma = zeros(gxvnum,gxvnum,max(max_ptildel));
        
        % if feedbacks are in, OLS coefficients related to the lagged
        % changes of feedbacks
        if not(isempty(exog_du))
            % compute the maximum number of feedbacks across global
            % variables
            maxfeednum = 0; feedbacks_flagmat_redux = [];
            for j=1:cols(feedbacks_flagmat)
                pointer = sum(feedbacks_flagmat(:,j));
                if pointer > 0
                    maxfeednum = maxfeednum + 1;
                    feedbacks_flagmat_redux = [feedbacks_flagmat_redux feedbacks_flagmat(:,j)]; %#ok
                end
            end
            
            Gammatilde = zeros(gxvnum,maxfeednum,max(max_qtildel));
        end
        
        detind = 1;
        for g = 1:gxvnum
            
            deplab_tmp = sprintf('d%s',gxvnames{g});
            deplab_tmp = {deplab_tmp};
            deplab = deplab_tmp;
            
            detlab = ['Intercept']; %#ok
            
            if not(estcase_du == 2)
            a0_du(g) = olsest.(gxvnames{g})(detind);
            end
            
            ptildel = max_ptildel(g);
            qtildel = max_qtildel(g);
            
           
            xlab = [];
            k=0;
            for i=1:ptildel
                for j = 1:gxvnum
                    Gamma(g,j,i) = olsest.(gxvnames{g})(detind+k+1);
                    k=k+1;
                    
                    xlab_tmp = sprintf('d%s_%d',gxvnames{j},i);
                    xlab_tmp = {xlab_tmp};
                    xlab = [xlab xlab_tmp]; %#ok
                end
            end
            
            xtildelab = [];
            if not(isempty(exog_du.(gxvnames{g})))
                k = 1;
                for i=1:qtildel
                    aj=1; 
                    for j=1:cols(feedbacks_flagmat_redux)

                        if feedbacks_flagmat_redux(g,j) == 1
                        Gammatilde(g,j,i) = olsest.(gxvnames{g})(detind+(ptildel*gxvnum)+k);
                        k=k+1;
                        
                        xtildelab_tmp = sprintf('d%s論%d',exognames_du.(gxvnames{g}){aj},i);
                        xtildelab_tmp = {xtildelab_tmp};
                        xtildelab = [xtildelab xtildelab_tmp]; %#ok                        
                        
                        aj=aj+1;
                        end
                    end    
                end
            end
            
            ecmlab = [];
            if not(rank_du == 0)
                for i=1:rank_du
                    ecmlab_tmp = sprintf('ecm%d_1',i);
                    ecmlab_tmp = {ecmlab_tmp};
                    ecmlab = [ecmlab ecmlab_tmp]; %#ok
                end
            end
            
            hlab = [deplab detlab ecmlab xlab xtildelab];
            hres = [NaN est_det.(gxvnames{g}) est_ecm.(gxvnames{g})' est_oth.(gxvnames{g})'];
            
            tab = [tab; hlab num2cell(NaN(1,ncells-length(hlab))); num2cell(hres) num2cell(NaN(1,ncells-length(hres)))]; %#ok
        end
        

        % recover the VARX coefficients of the dominant unit model
        %**********************************************************

        % - the Theta coeffs
        if not(isempty(ecm_du))
            Pi = alpha_du*Beta'; % Cointegration matrix
            sumthetasminuseye = Pi;
        else
            sumthetasminuseye = zeros(gxvnum);
        end
            
        if ptilde==1
            Theta_du(:,:,ptilde) = sumthetasminuseye + eye(gxvnum);
        elseif ptilde==2
            Theta_du(:,:,ptilde) = -Gamma(:,:,ptilde-1);
            Theta_du(:,:,1) = sumthetasminuseye +Gamma(:,:,1) + eye(gxvnum);
        else
            Theta_du(:,:,ptilde) = -Gamma(:,:,ptilde-1);
            for j=1:ptilde-2
                jj=ptilde-j;
                Theta_du(:,:,jj) = -(Gamma(:,:,jj-1)-(Gamma(:,:,jj)));
            end
            Theta_du(:,:,1) = sumthetasminuseye + Gamma(:,:,1) + eye(gxvnum);
        end
        
        if ptilde < amaxlag
            for j=1:amaxlag-ptilde
                Theta_du(:,:,ptilde+j) = zeros(gxvnum,gxvnum);
            end
        end

        
        if isfeedback == 1
            % - the Lambda coefficients
   
        if qtilde==1
        % This case will never arise       
        elseif qtilde==2
            Lambda_du(:,:,qtilde) = -Gammatilde(:,:,qtilde-1);
            Lambda_du(:,:,1) = Gammatilde(:,:,1);
        else
            Lambda_du(:,:,qtilde) = -Gammatilde(:,:,qtilde-1);
            for j=1:qtilde-2
                jj=qtilde-j;
                Lambda_du(:,:,jj) = -(Gammatilde(:,:,jj-1)-(Gammatilde(:,:,jj)));       
            end
            Lambda_du(:,:,1) = Gammatilde(:,:,1);
        end
         
            if qtilde < amaxlag
                for j=1:amaxlag-qtilde
                    Lambda_du(:,:,qtilde+j) = zeros(gxvnum,maxfeednum);
                end
            end
        else
            Lambda_du = [];
        end
    end
    
    % add the F-stats
    tab = [tab;  num2cell(NaN(3,ncells))];
    tab = [tab; 'F-Statistics for the Serial Correlation Test of Residuals' num2cell(NaN(1,ncells-1))];
    tab = [tab;  num2cell(NaN(1,ncells))];
    hlab = [' ' 'Fcrit_0.05' gxvnames];
    tab = [tab; hlab num2cell(NaN(1,ncells-cols(hlab)))];
    flabel_tmp = sprintf('F(%d,%d)', Fcommon(1),Fcommon(2));
    flabel_tmp = {flabel_tmp};
    ftab = [flabel_tmp num2cell([Fcommon(3) Feqspec]) num2cell(NaN(1,ncells-1 -cols([Fcommon(3) Feqspec])))];
    tab =[tab; ftab num2cell(NaN(1,ncells-cols(ftab)))];
    
else
    a0_du = [];
    a1_du = [];
    Theta_du = [];
    Lambda_du = [];
    
    
    % print the results of the lag order selection
    ncells = 10;
    tab = ['Choice Criteria for Selecting the Order of the Augmented Regression(s) Together With Corresponding Residual Serial Correlation F-Statistics' num2cell(NaN(1,ncells-1))];
    tab = [tab; num2cell(NaN(1,ncells))];
    tab = [tab; {' ' 'ptilde' 'qtilde' 'AIC' 'SBC' 'logLik' ' ' ' ' 'Fcrit_0.05' 'Fstat'}];
    tab = [tab; vlabel num2cell(vres) num2cell(NaN(rows(vres),1)) flabel num2cell(fres)]; 
end