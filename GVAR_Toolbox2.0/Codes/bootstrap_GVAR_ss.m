function  [kpsup_90cv kpmnsq_90cv ny_90cv rny_90cv qlr_90cv mw_90cv apw_90cv rqlr_90cv rmw_90cv rapw_90cv...
           kpsup_95cv kpmnsq_95cv ny_95cv rny_95cv qlr_95cv mw_95cv apw_95cv rqlr_95cv rmw_95cv rapw_95cv ...
           kpsup_99cv kpmnsq_99cv ny_99cv rny_99cv qlr_99cv mw_99cv apw_99cv rqlr_99cv rmw_99cv rapw_99cv] = bootstrap_GVAR_ss(ss_B,...
           ccut,kpsup,kpmnsq,ny,rny,qlr,mw,apw,rqlr,rmw,rapw,zeta,cb_varcov,y,ynames,cnum,cnames,...
          delta_0,delta_1,C,H0,dvtype,ntypes,wmatrices,maxlag,mlag,varxlag,estcase,rank,...
          dv, vnames, vnames_long, fvnames, gvnames, gvnames_long, vnum, gvnum, gxidx, gnvnum,gnvnames, ...
          dvflag, fvflag, gvflag,ydate,ylabelseq_estimation,yearseq_estimation,nyears_estimation,nobs,overid_flag,beta_r,shuffleflag)

%**************************************************************************
% PURPOSE: Performs the bootstrap of the GVAR model for computing the 
% critical values of the structural stability tests
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
    
    % do choleski
    P = chol(cb_varcov);
    Acap =P';
    invAcap = eye(rows(Acap))/Acap;
    
    eta = invAcap*zeta_c;
    etav = vec(eta);
    
elseif shuffleflag==1
    etav =zeta_c;
end

% Nonparametric approach
rand('state',8675432); %#ok


% Bootstrap GVAR
%***************

invH0 = eye(rows(H0))/H0;

kpsupres = struct(kpsup);
kpmnsqres = struct(kpmnsq);
nyres = struct(ny);
rnyres = struct(rny);
qlrres = struct(qlr);
mwres = struct(mw);
apwres = struct(apw);
rqlrres = struct(rqlr);
rmwres = struct(rmw);
rapwres = struct(rapw);

allvnames = [vnames gvnames];

packstep = 20;
packseq = 1:packstep:ss_B; % sequence for packing memory of MatLab, every 20 replications 
% the pack command is performed, see help
packseqidx = 1;

for b=1:ss_B
    
    y_b=zeros(rows(y),cols(y)); % preallocating matrices of bootstrapped endogenous variables
    zeta_bs = zeros(rows(y),cols(y));
    
    
    if b == packseq(packseqidx)
        disp('- Packing stored objects in workspace');
        pack storevars
        if not(b>ss_B-packstep)
            packseqidx = packseqidx +1;
        end
    end    
    
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

    % Creating y_t(b) series
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
    
    if not(gvnum==0)
        % need to disentangle whether each global variable is endogenous or exogenous to the
        % GVAR
        gxv_b = [];
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
    else
        gv_b = [];
    end
    
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
        for k=1:ntypes
            wmattype = sprintf('wmat%d',k);
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

    % create country models
    %**********************
    [endog_b endoglist junk exog_b junk junk junk] = create_countrymodels(cnames,cnum, ...
        vnames,vnames_long,fvnames,vnum,gvnames,gvnames_long,gvnum,gnvnum,gnvnames,dv_b,fv_b,gv_b,dvflag,fvflag,gvflag); %#ok
    
    
    % Estimate VECMX* (exactly identified estimation)
    %**********************************************
    for n = 1: cnum
        [beta_b_tmp.(cnames{n}) alpha_b.(cnames{n}) Psi_b.(cnames{n}) epsilon_b.(cnames{n})...
            Omega_b.(cnames{n}) ecm_b.(cnames{n})] = mlcoint(maxlag,rank.(cnames{n}), ...
            estcase.(cnames{n}),[],endog_b.(cnames{n}),varxlag.(cnames{n})(1), exog_b.(cnames{n}),varxlag.(cnames{n})(2));
    end
    
    if overid_flag == 1
        % Do overidentifying restriction estimation
        %******************************************
        
        for n=1:cnum
            if isfield(beta_r,cnames(n))
                
                [alpha_b.(cnames{n}) Psi_b.(cnames{n}) epsilon_b.(cnames{n}) ...
                    Omega_b.(cnames{n}) ecm_b.(cnames{n})] = mlcoint_r(beta_r.(cnames{n}),...
                    maxlag,estcase.(cnames{n}),endog_b.(cnames{n}),varxlag.(cnames{n})(1),...
                    exog_b.(cnames{n}),varxlag.(cnames{n})(2));
                
                beta_b_tmp.(cnames{n}) = beta_r.(cnames{n});
            end
        end
    end
    


    % Structural stability tests
    %***************************
    [kpsup_b kpmnsq_b ny_b rny_b qlr_b mw_b apw_b rqlr_b rmw_b rapw_b] = structural_stability_tests(cnum,cnames,endoglist,endog_b,exog_b,varxlag,estcase,maxlag,ecm_b,ccut);
    
    
    % store results of structural stability tests
    for i=1:length(allvnames)
        if isfield(kpsup,allvnames{i})
            clist = fieldnames(kpsup.(allvnames{i}));
            for n=1:length(clist);
                if b==1
                    kpsupres.(allvnames{i}).(clist{n}) = kpsup_b.(allvnames{i}).(clist{n});
                    kpmnsqres.(allvnames{i}).(clist{n}) = kpmnsq_b.(allvnames{i}).(clist{n});
                    nyres.(allvnames{i}).(clist{n}) = ny_b.(allvnames{i}).(clist{n});
                    rnyres.(allvnames{i}).(clist{n}) = rny_b.(allvnames{i}).(clist{n});
                    qlrres.(allvnames{i}).(clist{n}) = qlr_b.(allvnames{i}).(clist{n});
                    mwres.(allvnames{i}).(clist{n}) = mw_b.(allvnames{i}).(clist{n});
                    apwres.(allvnames{i}).(clist{n}) = apw_b.(allvnames{i}).(clist{n});
                    rqlrres.(allvnames{i}).(clist{n}) = rqlr_b.(allvnames{i}).(clist{n});
                    rmwres.(allvnames{i}).(clist{n}) = rmw_b.(allvnames{i}).(clist{n});
                    rapwres.(allvnames{i}).(clist{n}) = rapw_b.(allvnames{i}).(clist{n});
                else
                    kpsupres.(allvnames{i}).(clist{n}) = [kpsupres.(allvnames{i}).(clist{n}) kpsup_b.(allvnames{i}).(clist{n})];
                    kpmnsqres.(allvnames{i}).(clist{n}) = [kpmnsqres.(allvnames{i}).(clist{n}) kpmnsq_b.(allvnames{i}).(clist{n})];
                    nyres.(allvnames{i}).(clist{n}) = [nyres.(allvnames{i}).(clist{n}) ny_b.(allvnames{i}).(clist{n})];
                    rnyres.(allvnames{i}).(clist{n}) = [rnyres.(allvnames{i}).(clist{n}) rny_b.(allvnames{i}).(clist{n})];
                    qlrres.(allvnames{i}).(clist{n}) = [qlrres.(allvnames{i}).(clist{n}) qlr_b.(allvnames{i}).(clist{n})];
                    mwres.(allvnames{i}).(clist{n}) = [mwres.(allvnames{i}).(clist{n}) mw_b.(allvnames{i}).(clist{n})];
                    apwres.(allvnames{i}).(clist{n}) = [apwres.(allvnames{i}).(clist{n}) apw_b.(allvnames{i}).(clist{n})];
                    rqlrres.(allvnames{i}).(clist{n}) = [rqlrres.(allvnames{i}).(clist{n}) rqlr_b.(allvnames{i}).(clist{n})];
                    rmwres.(allvnames{i}).(clist{n}) = [rmwres.(allvnames{i}).(clist{n}) rmw_b.(allvnames{i}).(clist{n})];
                    rapwres.(allvnames{i}).(clist{n}) = [rapwres.(allvnames{i}).(clist{n}) rapw_b.(allvnames{i}).(clist{n})];
                end
            end
        end
    end
    clc
end


% Obtain bootstrap critical values of structural stability tests
%***************************************************************

kpsup_90cv = struct(kpsup); kpsup_95cv = struct(kpsup); kpsup_99cv = struct(kpsup);
kpmnsq_90cv = struct(kpmnsq); kpmnsq_95cv = struct(kpmnsq); kpmnsq_99cv = struct(kpmnsq);
ny_90cv = struct(ny); ny_95cv = struct(ny); ny_99cv = struct(ny);
rny_90cv = struct(rny); rny_95cv = struct(rny); rny_99cv = struct(rny);
qlr_90cv = struct(qlr); qlr_95cv = struct(qlr); qlr_99cv = struct(qlr);
mw_90cv = struct(mw); mw_95cv = struct(mw); mw_99cv = struct(mw);
apw_90cv = struct(apw); apw_95cv = struct(apw); apw_99cv = struct(apw);
rqlr_90cv = struct(rqlr); rqlr_95cv = struct(rqlr); rqlr_99cv = struct(rqlr);
rmw_90cv = struct(rmw); rmw_95cv = struct(rmw); rmw_99cv = struct(rmw);
rapw_90cv = struct(rapw); rapw_95cv = struct(rapw); rapw_99cv = struct(rapw);

for i=1:length(allvnames)
    if isfield(kpsup,allvnames{i})
        clist = fieldnames(kpsup.(allvnames{i}));
        for n=1:length(clist);
            % 90% significance level
            kpsup_90cv.(allvnames{i}).(clist{n}) = quantile(kpsupres.(allvnames{i}).(clist{n}),0.90);
            kpmnsq_90cv.(allvnames{i}).(clist{n}) = quantile(kpmnsqres.(allvnames{i}).(clist{n}),0.90);
            ny_90cv.(allvnames{i}).(clist{n}) = quantile(nyres.(allvnames{i}).(clist{n}),0.90);
            rny_90cv.(allvnames{i}).(clist{n}) = quantile(rnyres.(allvnames{i}).(clist{n}),0.90);
            qlr_90cv.(allvnames{i}).(clist{n}) = quantile(qlrres.(allvnames{i}).(clist{n}),0.90);
            mw_90cv.(allvnames{i}).(clist{n}) = quantile(mwres.(allvnames{i}).(clist{n}),0.90);
            apw_90cv.(allvnames{i}).(clist{n}) = quantile(apwres.(allvnames{i}).(clist{n}),0.90);
            rqlr_90cv.(allvnames{i}).(clist{n}) = quantile(rqlrres.(allvnames{i}).(clist{n}),0.90);
            rmw_90cv.(allvnames{i}).(clist{n}) = quantile(rmwres.(allvnames{i}).(clist{n}),0.90);
            rapw_90cv.(allvnames{i}).(clist{n}) = quantile(rapwres.(allvnames{i}).(clist{n}),0.90);

            % 95% significance level
            kpsup_95cv.(allvnames{i}).(clist{n}) = quantile(kpsupres.(allvnames{i}).(clist{n}),0.95);
            kpmnsq_95cv.(allvnames{i}).(clist{n}) = quantile(kpmnsqres.(allvnames{i}).(clist{n}),0.95);
            ny_95cv.(allvnames{i}).(clist{n}) = quantile(nyres.(allvnames{i}).(clist{n}),0.95);
            rny_95cv.(allvnames{i}).(clist{n}) = quantile(rnyres.(allvnames{i}).(clist{n}),0.95);
            qlr_95cv.(allvnames{i}).(clist{n}) = quantile(qlrres.(allvnames{i}).(clist{n}),0.95);
            mw_95cv.(allvnames{i}).(clist{n}) = quantile(mwres.(allvnames{i}).(clist{n}),0.95);
            apw_95cv.(allvnames{i}).(clist{n}) = quantile(apwres.(allvnames{i}).(clist{n}),0.95);
            rqlr_95cv.(allvnames{i}).(clist{n}) = quantile(rqlrres.(allvnames{i}).(clist{n}),0.95);
            rmw_95cv.(allvnames{i}).(clist{n}) = quantile(rmwres.(allvnames{i}).(clist{n}),0.95);
            rapw_95cv.(allvnames{i}).(clist{n}) = quantile(rapwres.(allvnames{i}).(clist{n}),0.95);

            % 99% significance level
            kpsup_99cv.(allvnames{i}).(clist{n}) = quantile(kpsupres.(allvnames{i}).(clist{n}),0.99);
            kpmnsq_99cv.(allvnames{i}).(clist{n}) = quantile(kpmnsqres.(allvnames{i}).(clist{n}),0.99);
            ny_99cv.(allvnames{i}).(clist{n}) = quantile(nyres.(allvnames{i}).(clist{n}),0.99);
            rny_99cv.(allvnames{i}).(clist{n}) = quantile(rnyres.(allvnames{i}).(clist{n}),0.99);
            qlr_99cv.(allvnames{i}).(clist{n}) = quantile(qlrres.(allvnames{i}).(clist{n}),0.99);
            mw_99cv.(allvnames{i}).(clist{n}) = quantile(mwres.(allvnames{i}).(clist{n}),0.99);
            apw_99cv.(allvnames{i}).(clist{n}) = quantile(apwres.(allvnames{i}).(clist{n}),0.99);
            rqlr_99cv.(allvnames{i}).(clist{n}) = quantile(rqlrres.(allvnames{i}).(clist{n}),0.99);
            rmw_99cv.(allvnames{i}).(clist{n}) = quantile(rmwres.(allvnames{i}).(clist{n}),0.99);
            rapw_99cv.(allvnames{i}).(clist{n}) = quantile(rapwres.(allvnames{i}).(clist{n}),0.99);
        end
    end
end

