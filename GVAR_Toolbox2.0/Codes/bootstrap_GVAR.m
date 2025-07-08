
function  [median_PP lbound_PP ubound_PP median_IRF lbound_IRF ubound_IRF median_FEVD lbound_FEVD ubound_FEVD ...
          ind_sim overid_LR_95cv overid_LR_99cv] = bootstrap_GVAR(B,N,zeta,cb_meth,lambda_star,use_shrinkedvcv,cb_varcov_dg, ...
          cb_country_exc,cnames_s_y,y,ynames,cnamesy,cnamesy_long,cnum,cnames,cweights,rnames,regions,rweights,...
          delta_0,delta_1,C,H0,K,dvtype,ntypes,wmatrices,maxlag,amaxlag,mlag,varxlag,estcase,rank,beta_r,...
          dv,vnames,vnames_long,fvnames,gvnames,gvnames_long,vnum,fvnum,gvnum,gnvnum,gnvnames,gxvnames,gxidx,gxcount,gxpointer,maxlag_du,varxlag_du,rank_du,estcase_du,...
          isfeedback,feedbacks_flagmat,fweights,exog_du,exognames_du,esttype_du2,estcase_du2,psc_du2,ptildel_found,qtildel_found,Wtilde,dvflag,fvflag,gvflag,allvnames,shockmtrx,...
          ydate,ylabelseq_estimation,yearseq_estimation,nyears_estimation,nobs,wmat_sol,redoflag,sgirfflag,firstcountries,newordervars,overid_flag,shuffleflag,yorder,nyorder)

%**************************************************************************
% PURPOSE: Performs the bootstrap of the GVAR model            
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


% Recenter residuals zeta
%*********************
zeta_c = zeros(rows(zeta),cols(zeta));
for i=1:rows(zeta)
    zeta_c(i,:) = zeta(i,:) - mean(zeta(i,:));
end

% Generate draws
%****************

if shuffleflag==0
    
    % the varcov matrix is positive definite, do Choleski
    P = chol(cb_varcov_dg);
    Acap =P';
    invAcap = eye(rows(Acap))/Acap;
    
    eta = invAcap*zeta_c;
    etav = vec(eta);
    
elseif shuffleflag==1
    etav =zeta_c;
end



% Nonparametric approach
rand('state',8675432);  %#ok

% Bootstrap GVAR
%***************

invH0 = eye(rows(H0))/H0;

PPres = [];
IRFres = [];
FEVDres = [];
if overid_flag == 1
    orcnames = fieldnames(beta_r);
    for n=1:length(orcnames)
        overid_LRres.(orcnames{n}) = [];
    end
end

packstep = 20;
packseq = 1:packstep:B; % sequence for packing memory of MatLab, every 20 replications
% the pack command is performed, see help
packseqidx = 1;

discdraws = 0; % counter for number of discarded draws
maxB = 2*B; % cap on total replications
reachedmax = 0;

for b=1:B
    
    y_b=zeros(rows(y),cols(y)); % preallocating matrices of bootstrapped endogenous variables
    zeta_bs = zeros(rows(y),cols(y));
    
    if b == packseq(packseqidx)
        disp('- Packing stored objects in workspace');
        pack storevars
        if not(b>B-packstep)
            packseqidx = packseqidx +1;
        end
    end
    
    k=1;
    % if model unstable, retry boostrap several times
    while k<2 % this conditional loop is used for checking stability of the bootstrapped GVAR model
        
        if shuffleflag==0
            ide = ceil(rand(rows(etav),1)*rows(etav));
            etav_b  = etav(ide,:);
            eta_b = reshape(etav_b,cols(eta),rows(eta))';
            zeta_bs(:,mlag+1:end) = Acap*eta_b;
            % first mlag observations are empty (as bootstrapping GVAR(mlag))
            % by premultiplying by Acap (the P) we are putting back the correlation
            % structure of the residuals; zeta_bs are the bootstrapped zeta
        elseif shuffleflag==1
            ide = ceil(rand(cols(etav),1)*cols(etav));
            etav_b  = etav(:,ide);
            zeta_bs(:,mlag+1:end) =etav_b;
        end
        
        eps = invH0*zeta_bs;
        
        fprintf(' Bootstrap #:%-2d \r',b);
        
        disp('- Creating the bootstrap series');
        % Creating x_t(b) series
        %***********************
        y_b(:,1:mlag) = y(:,1:mlag); % maintain first mlag observations
        for t=mlag+1:cols(y)
            acc=0;
            for j=1:mlag
                acc = acc + C(:,:,j)*y_b(:,t-j);
            end
            y_b(:,t) = delta_0 +delta_1*(t-1) + acc + eps(:,t);
        end
        
        % Retrieve domestic variables
        %****************************
        for j=1:vnum
            varlist = fieldnames(dv.(vnames{j}));
            varn = 1;
            for i=1:rows(y)
                if strcmp(ynames{i},vnames{j})
                    dv_b.(vnames{j}).(varlist{varn}) = y_b(i,:)';
                    varn=varn+1;
                end
            end
        end
        
        gv_b = [];
        gxv_b = [];
        if not(gvnum==0)
            % need to disentangle whether each global variable is endogenous or exogenous to the
            % GVAR
            for g=1:gvnum
                for i=1:rows(y)
                    if strcmp(ynames{i},gvnames{g});
                        gv_b.(gvnames{g}) = y_b(i,:)';
                    end
                end
                if not(isempty(gxidx)) && gxidx(g) == 1
                    gxv_b = [gxv_b gv_b.(gvnames{g})]; %#ok
                end
            end
        end
        gxvnum = cols(gxv_b);
        
        % Create foreign-specific variables
        %**********************************
        for t=1:nyears_estimation
            
            % generate the structure of domestic variables
            obs_idx = [];
            for j=1:nobs
                if  ydate(j) == yearseq_estimation(t)
                    obs_idx = [obs_idx; j]; %#ok
                end
            end
            for i=1:vnum
                clist = fieldnames(dv_b.(vnames{i}));
                for n=1:length(clist)
                    dv_b_t.(ylabelseq_estimation{t}).(vnames{i}).(clist{n}) = dv_b.(vnames{i}).(clist{n})(obs_idx);
                end
            end
            
            % generate the structure of weight matrices
            for kk=1:ntypes
                wmattype = sprintf('wmat%d',kk);
                wm_t.(ylabelseq_estimation{t}).(wmattype) = wmatrices.(wmattype)(:,:,t);
            end
            
            % generate the structure of foreign variables
            fv_b_t.(ylabelseq_estimation{t}) = create_foreignvariables(cnum, cnames, vnum, vnames, dv_b_t.(ylabelseq_estimation{t}),dvtype, wm_t.(ylabelseq_estimation{t}));
        end
        
        % now assemble the series of foreign-specific variables
        for t=1:nyears_estimation
            if t==1
                fv_b = struct(fv_b_t.(ylabelseq_estimation{t}));
            else
                % reconstruct the series:
                for i=1:vnum
                    clist = fieldnames(fv_b_t.(ylabelseq_estimation{t}).(vnames{i}));
                    for n=1:length(clist)
                        fv_b.(vnames{i}).(clist{n}) = [fv_b.(vnames{i}).(clist{n}); fv_b_t.(ylabelseq_estimation{t}).(vnames{i}).(clist{n})];
                    end
                end
            end
        end
        
        
        % cleaning variables
        clear fv_b_t
        
        
        if not(isempty(gxv_b)) % Dominant unit model
            if gxvnum > 1 % multivariate model
                
                % VECM estimation
                % skip if there are no coefficients to estimate
                if (rank_du == 0) && (estcase_du == 2) && (varxlag_du(1)  == 1)
                    % there are no coefficients to estimate in the VEC model        
                    ecm_du = [];
                    alpha_du = [];
                    beta_du = [];
                    beta_b.du_model = [];
                    beta_norm_b.du_model = [];       
                else
                    [beta_du alpha_du Psi_du junk junk ecm_du junk junk junk junk junk junk ...
                        junk junk junk junk junk] = mlcoint(maxlag_du,...
                        rank_du,estcase_du,[],gxv_b,varxlag_du(1)); %#ok
                    
                    % normalize the beta vector
                    if estcase_du == 4 || estcase_du == 2 % exclude trend or intercept from beta vector
                        invb = eye(rows(beta_du(2:rank_du+1,1:rank_du)))/beta_du(2:rank_du+1,1:rank_du);
                    elseif estcase_du == 3
                        invb = eye(rows(beta_du(1:rank_du,1:rank_du)))/beta_du(1:rank_du,1:rank_du);
                    end
                    beta_norm_du = beta_du*invb;
                    
                    % store beta_du and beta_norm_du in the beta and beta_norm
                    % structures (which contain the beta vectors of country models)
                    beta_b.du_model = beta_du;
                    beta_norm_b.du_model = beta_norm_du;
                end
            else
                ecm_du = [];
                alpha_du = [];
                beta_du = [];
                estcase_du = [];
                beta_norm_b.du_model = [];
            end
            
            
            if isfeedback == 1
                % build feedback variables
                for g = 1:gxvnum
                    
                    exog_du.(gxvnames{g}) = [];
                    
                    % identify which feedback variables will enter in equation of
                    % global variable g
                    fvdunames.(gxvnames{g}) = {};
                    for kk=1:length(vnames)
                        if feedbacks_flagmat(g,kk) == 1
                            fvdunames.(gxvnames{g}) = [fvdunames.(gxvnames{g}) vnames(kk)];
                        end
                    end
                    
                    % create the feedback variables
                    fvlist = fvdunames.(gxvnames{g});
                    for j = 1 : length(fvlist)
                        clist = fieldnames(dv.(fvlist{j}));
                        
                        sumx = 0;
                        
                        for i = 1:length(clist)
                            sumx = sumx + fweights.(fvlist{j}).(clist{i})*dv_b.(fvlist{j}).(clist{i});
                        end
                        
                        fvdu_b.(fvlist{j}) = sumx;
                        exog_du.(gxvnames{g}) = [exog_du.(gxvnames{g}) fvdu_b.(fvlist{j})];
                    end
                    
                    exognames_du.(gxvnames{g}) = fvdunames.(gxvnames{g});
                end
            end
            
            
            if isfeedback == 0 && gxvnum > 1
                % no augmented regression, recover the VAR coefficients of the
                % (multivariate) dominant unit model from the VECM estimates
                [a0_du_b a1_du_b Theta_du_b] =  vec2var_du(amaxlag,gxvnum,varxlag_du,alpha_du,beta_du,Psi_du,estcase_du);
                maxfeednum = length(vnames);
                Lambda_du_b = zeros(gxvnum,maxfeednum,amaxlag);
            else
                % Estimation of the augmented regression
                [junk junk junk a0_du_b a1_du_b Theta_du_b Lambda_du_b junk] = augmentedregression(isfeedback,0,[],...
                    ptildel_found,qtildel_found,maxlag,gxv_b,gxvnames,exog_du,exognames_du,feedbacks_flagmat,ecm_du,alpha_du,beta_du,estcase_du,esttype_du2,estcase_du2,psc_du2);  %#ok
                         
            end          
        else
            a0_du_b = [];
            a1_du_b = [];
            Theta_du_b = [];
            Lambda_du_b = [];
            beta_norm_b.du_model = [];
        end
        
        
        
        % create country models
        %**********************
        [endog_b endoglist junk exog_b exoglist exognx_b exognxlist] = create_countrymodels(cnames,cnum, ...
            vnames,vnames_long,fvnames,vnum,gvnames,gvnames_long,gvnum,gnvnum,gnvnames,dv_b,fv_b,gv_b,dvflag,fvflag,gvflag); %#ok
        
        disp('- Estimating the VECMX* individual models ');
        % Estimate VECMX* (exactly identified estimation)
        %**********************************************
        for n = 1: cnum
            [beta_b_tmp.(cnames{n}) alpha_b.(cnames{n}) Psi_b.(cnames{n}) epsilon_b.(cnames{n})...
                Omega_b.(cnames{n}) ecm_b.(cnames{n}) std_b.(cnames{n}) logl_b.(cnames{n})] = mlcoint(maxlag,rank.(cnames{n}), ...
                estcase.(cnames{n}),[],endog_b.(cnames{n}),varxlag.(cnames{n})(1), exog_b.(cnames{n}),varxlag.(cnames{n})(2));
        end
        
        if overid_flag == 1
            % Do overidentifying restriction estimation
            %******************************************
            
            for n=1:cnum
                
                if isfield(beta_r,cnames(n))
                    
                    [alpha_b.(cnames{n}) Psi_b.(cnames{n}) epsilon_b.(cnames{n}) ...
                        Omega_b.(cnames{n}) ecm_b.(cnames{n}) std_b.(cnames{n}) ...
                        logl_b_r.(cnames{n})] = mlcoint_r(beta_r.(cnames{n}),...
                        maxlag,estcase.(cnames{n}),endog_b.(cnames{n}),varxlag.(cnames{n})(1),...
                        exog_b.(cnames{n}),varxlag.(cnames{n})(2));
                    
                    beta_b_tmp.(cnames{n}) = beta_r.(cnames{n});
                    
                    overid_LR.(cnames{n}) = -2*(logl_b_r.(cnames{n})-logl_b.(cnames{n}));
                    
                    % store results
                    overid_LRres.(cnames{n}) = [overid_LRres.(cnames{n}) overid_LR.(cnames{n})];
                end
                
                
            end
        end
        
        % Creating link matrices
        %***********************
        % use wmat_sol, which has previously computed (both for fixed and
        % time-varying weights)
        [z_b x_b junk W_b] = create_linkmatrices(cnames,cnum,vnum,vnames,fvnum, fvnames, gnvnum,gnvnames,endog_b,endoglist,exognx_b,exognxlist,dvtype,wmat_sol); %#ok
        
        % create the y vector, which includes the x vector and the global exogenous variables
        y_b = [x_b; gxv_b'];
        
        if not(isempty(gxv_b))
            % retrieve also the z vector augmented by the global (exogenous)
            % variables, denote it zy
            zy_b = create_linkmatrices(cnames,cnum,vnum,vnames,fvnum, fvnames, gnvnum,gnvnames,endog_b,endoglist,exog_b,exoglist,dvtype,wmat_sol);
            
            % augment the W matrices accordingly
            for n=1:cnum
                Wy_b.(cnames{n}) = [W_b.(cnames{n}) zeros(rows(W_b.(cnames{n})),gxvnum); zeros(gxvnum,K) eye(gxvnum)];
            end
            
            % create also an augmented W matrix for the dominant unit model
            Wy_b.du_model = [zeros(gxvnum,K) eye(gxvnum)];
        else
            zy_b = z_b;
            Wy_b = W_b;
        end
        
        % Retrieve VARX models parameters from VECMX estimates
        %****************************************************
        [a0_b a1_b Theta_b Lambda0_tmp_b Lambda_tmp_b] =  vecx2varx(maxlag,cnum,cnames,zy_b,endoglist,varxlag,alpha_b,beta_b_tmp,Psi_b,estcase);
        
        lsize = size(Lambda_tmp_b.(cnames{1}));
        if length(lsize) == 3
            nlag = lsize(3);
        elseif length(lsize) == 2
            nlag = 1;
        end
        
        
        for n=1:cnum
            Lambda0_b.(cnames{n}) = [];
            Lambda_b.(cnames{n}) = [];
            
            Gamma0_b.(cnames{n}) = zeros(length(endoglist.(cnames{n})),gxvnum);
            Gamma_b.(cnames{n}) = zeros(length(endoglist.(cnames{n})),gxvnum,nlag);
            
            for j=1:length(exoglist.(cnames{n}))
                flag = 0;
                for i=1:length(gxvnames)
                    if strcmp(exoglist.(cnames{n})(j),gxvnames{i})
                        flag = 1;
                        Gamma0_b.(cnames{n})(:,i) = Lambda0_tmp_b.(cnames{n})(:,j);
                        Gamma_b.(cnames{n})(:,i,:) = Lambda_tmp_b.(cnames{n})(:,j,:);
                    end
                end
                if flag == 0;
                    Lambda0_b.(cnames{n}) = [Lambda0_b.(cnames{n}) Lambda0_tmp_b.(cnames{n})(:,j)];
                    Lambda_b.(cnames{n}) = [Lambda_b.(cnames{n}) Lambda_tmp_b.(cnames{n})(:,j,:)];
                end
            end
        end
        
        
        % compute beta_vectors
        %*********************
        for n=1:cnum
            beta_b.(cnames{n}) = beta_b_tmp.(cnames{n});
            
            if overid_flag == 1 && isfield(beta_r,cnames(n)) % then don't do any normalization
                beta_norm_b.(cnames{n}) = beta_b.(cnames{n});
            else
                % normalize beta's
                rk = rank.(cnames{n});
                b_b = beta_b.(cnames{n});
                
                if estcase.(cnames{n}) == 4 || estcase.(cnames{n}) == 2% exclude trend / intercept from beta
                    invb = eye(rows(b_b(2:rk+1,1:rk)))/b_b(2:rk+1,1:rk);
                elseif estcase.(cnames{n}) == 3 
                    invb = eye(rows(b_b(1:rk,1:rk)))/b_b(1:rk,1:rk);
                end
                beta_norm_b.(cnames{n}) = b_b*invb;
            end
        end
        
        
        disp('- Solving the GVAR');
    
        % Solve GVAR
        %***********
        [junk junk H0_b C_b zeta_b Sigma_b] = solve_GVAR(maxlag,amaxlag,...
            cnum,cnames,W_b,a0_b,a1_b,Theta_b,Lambda0_b,Lambda_b,Gamma0_b,Gamma_b,x_b,gxv_b,Wtilde,a0_du_b,a1_du_b,...
            Theta_du_b,Lambda_du_b); %#ok
        % note: in the code above, we are using the W matrices (and not Wy, which are augmented for global exog variables)

        
        % calculate companion form
        %*************************
        Ky = K +gxvnum;
        
        Cbarup_b = [];
        for j=1:mlag
            Cbarup_b = [Cbarup_b C_b(:,:,j)]; %#ok
        end
        Cbardown_b = [eye(Ky*(mlag-1)) zeros(Ky*(mlag-1),Ky)];
        Cbar_b = [Cbarup_b; Cbardown_b];
        
        % calculating eigenvalues of the GVAR model
        %*******************************************
        eigens_b = eig(Cbar_b);
        
        tol = 0.0001;  % tolerance parameter
        if max(abs(eigens_b)) > 1+tol  % modified: before was max(eigens_b) > 1+tol. A.G. 2012
            disp('- GVAR is unstable: Another bootstrap replication will be performed');
            disp(max(abs(eigens_b)))
            discdraws = discdraws +1;
        else
            disp('- GVAR is stable (i.e. all eigenvalues are <=1)');
            k=2;
        end
        
        if b+discdraws == maxB
            reachedmax = 1; % the bootstrap cap is reached
            break
        end
    end
    
    if reachedmax == 1
        break
    end
    
    
    

    
    % compute varcov matrix
    %**********************
    % - do transformation
    cb_varcov_tx = transform_varcov(cb_meth,cb_country_exc,Sigma_b,cnames_s_y);
    if use_shrinkedvcv == 1
        % - do shrinkage
        cb_varcov = ShrinkageCorrLstar(cb_varcov_tx,rows(zeta),lambda_star);
    elseif use_shrinkedvcv == 0
        cb_varcov = cb_varcov_tx;
    end
    
    if sgirfflag == 1 || sgirfflag == 2
        % Reorder GVAR for Structural GIRF
        %*********************************
        
        [cb_varcov_s Sigma_zeta0_b H0_b_s_t C_b_s Wy_b_s beta_b_s beta_norm_b_s ncnames ncnames_long ...
            junk nynames junk ncnames_s_y] = reorder_GVAR(firstcountries,...
            newordervars,cnamesy,cnamesy_long,Ky,endoglist,...
            H0_b,C_b,mlag,ynames,cnames_s_y,Wy_b,beta_b,beta_norm_b,estcase,cb_varcov); %#ok
        
    else
        Sigma_zeta0_b = [];
        H0_b_s_t = [];
    end
    
    
    % calculating dynamic multipliers
    %********************************
    if sgirfflag == 1 || sgirfflag == 2
        PHI = dyn_multipliers(Ky,mlag,C_b,N);
        PHI_s = dyn_multipliers(Ky,mlag,C_b_s,N);
    else
        PHI = dyn_multipliers(Ky,mlag,C_b,N);
    end
    
    if redoflag == 0 % the first time you run the program, persistence profiles are computed
        disp('- Computing the PPs');
        % build PProfiles
        %****************
 
        if isempty(beta_norm_b.du_model)
            cnamesy_pp = cnamesy(not(strcmp(cnamesy,'du_model')));
        else
            cnamesy_pp = cnamesy;
        end
        
        if gxvnum > 0
            % augment beta_norm matrices with zeros for excluded global exogenous
            % variables, so that matrices are conformable in computing persistence
            % profiles
            for n=1:length(cnamesy)
                if not(strcmp(cnamesy(n),'du_model'))
                    % by construction the dominant unit model contains all global
                    % exogenous variables, so skip it
                    oldblock = beta_norm_b.(cnamesy{n})(end+1-gxcount.(cnamesy{n}):end,:);
                    newblock = zeros(gxvnum,rank.(cnamesy{n}));
                    kk = 1;
                    for k=1:gxvnum
                        if gxpointer.(cnamesy{n})(k) == 1
                            newblock(k,:) = oldblock(kk,:);
                            kk = kk + 1;
                        end
                    end
                    beta_norm_gx_b.(cnamesy{n}) = [beta_norm_b.(cnamesy{n})(1:end-gxcount.(cnamesy{n}),:); newblock];
                end
            end
            beta_norm_gx_b.du_model = beta_norm_b.du_model;
        else
            beta_norm_gx_b = beta_norm_b;
        end
        
        PP_temp = pprofile(PHI,cb_varcov,H0_b,Wy_b,beta_norm_gx_b,N,cnamesy_pp,estcase);
        
        % store PPprofiles
        PPres(:,b,:) = PP_temp; %#ok
    end
    
    
    disp('- Computing the IRFs & FEVDs');
    % IRFs and associated FEVDs
    %****************************
    ind_sim = 0;
    for j=1:cols(shockmtrx)
        for n=1:rows(shockmtrx)
            
            if not(shockmtrx(n,j) == 0) && not(isnan(shockmtrx(n,j)))
                ind_sim = ind_sim+1;
                
                eslct = zeros(length(ynames),1);   % selection vector
                
                shocksign = sign(shockmtrx(n,j));
                % if == 1: positive shock
                % if == -1: negative shock
                shocktype = abs(shockmtrx(n,j));
                % if == 1: country-specific shock
                % if == 2: region-specific shock
                % if == 3: global-shock
                if shocktype == 1  % country-specific shock
                    regtype = [];
                    countype = cnamesy(n);
                    vartype = allvnames(j);
                    for i=1:length(ynames)
                        if (strcmp(ynames{i},vartype) && strcmp(cnames_s_y{i},countype))
                            eslct(i) = 1;
                        end
                    end
                    
                elseif shocktype == 2
                    country = cnames(n);
                    for r=1:length(rnames)
                        for jj=1:length(regions.(rnames{r}))
                            if strcmp(country,regions.(rnames{r})(jj))
                                regtype = rnames(r);
                            end
                        end
                    end
                    
                    
                    vartype = allvnames(j);
                    for k = 1:length(regions.(regtype{1}))
                        for i=1:length(ynames)
                            if (strcmp(ynames{i},vartype) && strcmp(cnames_s_y{i},regions.(regtype{1}){k}))
                                eslct(i) = rweights.(vartype{1}).(regtype{1}).(regions.(regtype{1}){k});   % note: vartype{1} is only a trick
                            end
                        end
                    end
                    
                    
                elseif shocktype == 3  % global shock
                    
                    regtype=[];
                    vartype = allvnames(j);
                    countrylist = fieldnames(cweights.(allvnames{j}));
                    k=0;
                    for i=1:length(ynames)
                        if strcmp(ynames{i},vartype)
                            k=k+1;
                            eslct(i) = cweights.(vartype{1}).(countrylist{k});
                            % note: typing vartype{1} is only a trick to pass from string to char
                        end
                    end
                end
                
                if shocksign == -1
                    eslct = -eslct;
                end
                
                
                % build IRFs and FEVDs
                %**********************
                if sgirfflag == 1 || sgirfflag == 2 
                    
                    % update the shock selector given the new reordering of
                    % variables (amended 29 September 2014, A.G.)
                    oldpos = yorder(not(eslct==0));
                    for i=1:length(oldpos)
                        newpos = find(nyorder==oldpos(i));
                        eslct(oldpos(i)) = 0;
                        if shocksign == 1
                            eslct(newpos) = 1;
                        elseif shocksign == -1
                            eslct(newpos) = -1;
                        end
                    end
                    
                    % use reordered matrices PHI_s, cb_varcov_s, H0_b_s_t
                    IRFMAT = irf(Ky,N,PHI_s,cb_varcov_s,H0_b_s_t,eslct,sgirfflag,Sigma_zeta0_b);
                    FEVDMAT = fevd(Ky,N,PHI_s,cb_varcov_s,H0_b_s_t,eslct,sgirfflag,Sigma_zeta0_b);                 
                             
                else
                    IRFMAT = irf(Ky,N,PHI,cb_varcov,H0_b,eslct,sgirfflag,Sigma_zeta0_b);
                    FEVDMAT = fevd(Ky,N,PHI,cb_varcov,H0_b,eslct,sgirfflag,Sigma_zeta0_b);
                end
                
                % store IRFs and FEVDs
                IRFres(:,b,:,ind_sim) = IRFMAT; %#ok
                FEVDres(:,b,:,ind_sim) = FEVDMAT; %#ok
            end
        end
    end   
end

if reachedmax == 1
    actual = b+discdraws; % actual number of replications
    disp('Bootstrap completed (reaching the maximum number of bootstrap replications).');
else
    actual = B+discdraws; % actual number of replications
    disp('Bootstrap completed.');
end
fprintf('Number of actual draws:%-2d \r',actual);
disp(' ');

if overid_flag == 1
    % Obtain bootstrap critical values for overidentifying restrictions tests
    %**********************************************************************
    
    for n=1:cnum
        if isfield(beta_r,cnames(n))
            overid_LR_95cv.(cnames{n}) = quantile(overid_LRres.(cnames{n}),0.95);
            overid_LR_99cv.(cnames{n})= quantile(overid_LRres.(cnames{n}),0.99);
        else
            overid_LR_95cv.(cnames{n}) = NaN;
            overid_LR_99cv.(cnames{n})= NaN;
        end
    end
else
    overid_LR_95cv = [];
    overid_LR_99cv = [];
end

% Obtain 90% bootstrap confidence intervals and median estimates of PPs,
% generalized impulse response functions and FEVDs
%*************************************************************************

if redoflag == 0
    median_PP = NaN(rows(PPres),N+1);
    lbound_PP = NaN(rows(PPres),N+1);
    ubound_PP = NaN(rows(PPres),N+1);
    
    for i=1:rows(PPres)
        for t=1:N+1
            median_PP(i,t) = median(PPres(i,:,t));
            lbound_PP(i,t) = quantile(PPres(i,:,t),0.05);
            ubound_PP(i,t) = quantile(PPres(i,:,t),0.95);
        end
    end
    
else
    median_PP = [];
    lbound_PP = [];
    ubound_PP = [];
end

nsim = ind_sim; % number of shocks performed

median_IRF = NaN(rows(y),N+1,nsim);
lbound_IRF = NaN(rows(y),N+1,nsim);
ubound_IRF = NaN(rows(y),N+1,nsim);

median_FEVD = NaN(rows(y),N+1,nsim);
lbound_FEVD = NaN(rows(y),N+1,nsim);
ubound_FEVD = NaN(rows(y),N+1,nsim);

for h=1:nsim
    for i=1:rows(y)
        for t=1:N+1
            median_IRF(i,t,h) = median(IRFres(i,:,t,h));
            lbound_IRF(i,t,h) = quantile(IRFres(i,:,t,h),0.05);
            ubound_IRF(i,t,h) = quantile(IRFres(i,:,t,h),0.95);
            
            median_FEVD(i,t,h) = median(FEVDres(i,:,t,h));
            lbound_FEVD(i,t,h) = quantile(FEVDres(i,:,t,h),0.05);
            ubound_FEVD(i,t,h) = quantile(FEVDres(i,:,t,h),0.95);
        end
    end
end

