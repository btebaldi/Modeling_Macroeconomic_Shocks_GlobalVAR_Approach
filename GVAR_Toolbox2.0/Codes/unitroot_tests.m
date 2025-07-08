function [urt_adf_out urt_ws_out urt_adf urt_ws] = unitroot_tests(vnum,...
  fvnum,gvnum,vnames,gvnames,dv,fv,gv,maxlag_urt,info_urt_choice)

%**************************************************************************
% PURPOSE: performing unit root tests on domestic, foreign-specific and
% global variables (ADF and Weighted Symmetric Tests)
%--------------------------------------------------------------------------
% INPUT: 
% - vnum: number of domestic variables
% - fvnum: number of foreign-specific variables
% - gvnum: number of global variables
% - vnames: list of domestic variables
% - gvnames: list of global variables
% - dv: structure, contains data of domestic variables for each country
% - dv: structure, contains data of foreign-specific variables for each
% country
% - gv: structure, contains data of global variables for each country
% - maxlag_urt: maximum lag order of unit root regressions (the routine
% performs a set of regressions with increasing lag order up to maxlag_urt,
% then chooses the regression according to an information criterion 
% specified in info_urt_choice)
% - info_urt_choice: information criterion chosen for determining the lag 
% order of unit root tests regressions, it can be Akaike or Schwartz 
% Bayesian
% 
%--------------------------------------------------------------------------
% OUTPUT: 
% - urt_adf_out: structure, contains ADF stat, ADF 95% critical value and
% lag order employed in test, for each variable of each country
% - urt_ws_out: structure, contains WS stat, WS 95% critical value and
% lag order employed in test, for each variable of each country
% - urt_adf: structure, contains further info about ADF unit root tests 
% regressions, specifically Akaike, Schwartz Bayesian criteria and LogL 
% for each regression starting from lag order = 0 to lag order = urt_maxlag
% - urt_ws: structure, contains further info about ADF unit root tests 
% regressions, specifically Akaike, Schwartz Bayesian criteria and LogL 
% for each regression starting from lag order = 0 to lag order = urt_maxlag
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     


% asymptotic 5% significance critical values (see DdPS) 
adfcrit0 = -2.89; % - yes constant, no trend
adfcrit1 = -3.45; % - yes constant, yes trend
wscrit0 = -2.55; % - yes constant, no trend 
wscrit1 = -3.24; % - yes constant, yes trend

if strcmp(info_urt_choice,'aic')
    info_urt = 2;
elseif strcmp(info_urt_choice,'sbc')
    info_urt = 3;
end

warning off all


% deterministic components
c = 0; % intercept
ct = 1; % intercept plus trend

% domestic variables
disp('- Testing domestic variables')
for i=1:vnum
    clist = fieldnames(dv.(vnames{i}));
    for j=1:length(clist)
        % initialize output structures, in which the following items will be stored: adf/ws stat,
        % aic, sbc, logl        
        
        Levelname_t = sprintf('%s_t',vnames{i});
        urt_adf.(Levelname_t).(clist{j}) = zeros(maxlag_urt,4); % levels (with trend)
        urt_ws.(Levelname_t).(clist{j}) = zeros(maxlag_urt,4);
      
        Levelname_nt = sprintf('%s_nt',vnames{i});
        urt_adf.(Levelname_nt).(clist{j}) = zeros(maxlag_urt,4); % levels (no trend)
        urt_ws.(Levelname_nt).(clist{j}) = zeros(maxlag_urt,4);        
        
        Idiffname = sprintf('D%s',vnames{i}); 
        urt_adf.(Idiffname).(clist{j}) = zeros(maxlag_urt,4); % Ist diff
        urt_ws.(Idiffname).(clist{j}) = zeros(maxlag_urt,4);
        
        IIdiffname = sprintf('DD%s',vnames{i});
        urt_adf.(IIdiffname).(clist{j}) = zeros(maxlag_urt,4); % IInd diff
        urt_ws.(IIdiffname).(clist{j}) = zeros(maxlag_urt,4);
        
        for p=1:maxlag_urt
            % testing levels
            % with trend
            urt_adf.(Levelname_t).(clist{j})(p,:) = adf(dv.(vnames{i}).(clist{j}),ct,p);
            urt_ws.(Levelname_t).(clist{j})(p,:) = ws(dv.(vnames{i}).(clist{j}),ct,p);
            % without trend
            urt_adf.(Levelname_nt).(clist{j})(p,:) = adf(dv.(vnames{i}).(clist{j}),c,p);
            urt_ws.(Levelname_nt).(clist{j})(p,:) = ws(dv.(vnames{i}).(clist{j}),c,p);
            % testing Ist diffs  (just intercept included)
            urt_adf.(Idiffname).(clist{j})(p,:) = adf(diff(dv.(vnames{i}).(clist{j}),1),c,p);
            urt_ws.(Idiffname).(clist{j})(p,:) = ws(diff(dv.(vnames{i}).(clist{j}),1),c,p);
            % testing IInd diffs (just intercept included)
            urt_adf.(IIdiffname).(clist{j})(p,:) = adf(diff(dv.(vnames{i}).(clist{j}),2),c,p);
            urt_ws.(IIdiffname).(clist{j})(p,:) = ws(diff(dv.(vnames{i}).(clist{j}),2),c,p);
        end
        
        % selecting ADF results depending on an information criterion
        urt_adf_out.(vnames{i}).(clist{j}) = zeros(4,3);
        urt_ws_out.(vnames{i}).(clist{j}) = zeros(4,3);
        
        
        [junk posmin] = min(urt_adf.(Levelname_t).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(vnames{i}).(clist{j})(1,:) = [adfcrit1 urt_adf.(Levelname_t).(clist{j})(posmin,1) posmin];
        urt_ws_out.(vnames{i}).(clist{j})(1,:) = [wscrit1 urt_ws.(Levelname_t).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(Levelname_nt).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(vnames{i}).(clist{j})(2,:) = [adfcrit0 urt_adf.(Levelname_nt).(clist{j})(posmin,1) posmin];
        urt_ws_out.(vnames{i}).(clist{j})(2,:) = [wscrit0 urt_ws.(Levelname_nt).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(Idiffname).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(vnames{i}).(clist{j})(3,:) = [adfcrit0 urt_adf.(Idiffname).(clist{j})(posmin,1) posmin];
        urt_ws_out.(vnames{i}).(clist{j})(3,:) = [wscrit0 urt_ws.(Idiffname).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(IIdiffname).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(vnames{i}).(clist{j})(4,:) = [adfcrit0 urt_adf.(IIdiffname).(clist{j})(posmin,1) posmin]; 
        urt_ws_out.(vnames{i}).(clist{j})(4,:) = [wscrit0 urt_ws.(IIdiffname).(clist{j})(posmin,1) posmin];
               
    end
end



% foreign variables
disp('- Testing foreign-specific variables')
for i=1:fvnum
    clist = fieldnames(fv.(vnames{i}));
    for j=1:length(clist)
        % initialize output structures, in which the following items will be stored: adf/ws stat,
        % aic, sbc, logl        
        
        Levelname_t = sprintf('%ss_t',vnames{i});
        urt_adf.(Levelname_t).(clist{j}) = zeros(maxlag_urt,4); % levels (with trend)
        urt_ws.(Levelname_t).(clist{j}) = zeros(maxlag_urt,4);
        
        Levelname_nt = sprintf('%ss_nt',vnames{i});
        urt_adf.(Levelname_nt).(clist{j}) = zeros(maxlag_urt,4); % levels (without trend)
        urt_ws.(Levelname_nt).(clist{j}) = zeros(maxlag_urt,4);        
        
        Idiffname = sprintf('D%ss',vnames{i}); 
        urt_adf.(Idiffname).(clist{j}) = zeros(maxlag_urt,4); % Ist diff
        urt_ws.(Idiffname).(clist{j}) = zeros(maxlag_urt,4);
        
        IIdiffname = sprintf('DD%ss',vnames{i});
        urt_adf.(IIdiffname).(clist{j}) = zeros(maxlag_urt,4); % IInd diff
        urt_ws.(IIdiffname).(clist{j}) = zeros(maxlag_urt,4);
        
        for p=1:maxlag_urt
            % testing levels
            % with trend
            urt_adf.(Levelname_t).(clist{j})(p,:) = adf(fv.(vnames{i}).(clist{j}),ct,p);
            urt_ws.(Levelname_t).(clist{j})(p,:) = ws(fv.(vnames{i}).(clist{j}),ct,p);
            % without trend
            urt_adf.(Levelname_nt).(clist{j})(p,:) = adf(fv.(vnames{i}).(clist{j}),c,p);
            urt_ws.(Levelname_nt).(clist{j})(p,:) = ws(fv.(vnames{i}).(clist{j}),c,p);
            % testing Ist diffs  (just intercept included)
            urt_adf.(Idiffname).(clist{j})(p,:) = adf(diff(fv.(vnames{i}).(clist{j}),1),c,p);
            urt_ws.(Idiffname).(clist{j})(p,:) = ws(diff(fv.(vnames{i}).(clist{j}),1),c,p);
            % testing IInd diffs (just intercept included)
            urt_adf.(IIdiffname).(clist{j})(p,:) = adf(diff(fv.(vnames{i}).(clist{j}),2),c,p);
            urt_ws.(IIdiffname).(clist{j})(p,:) = ws(diff(fv.(vnames{i}).(clist{j}),2),c,p);
        end
        
        % selecting ADF results depending on an information criterion
        Levelname = sprintf('%ss',vnames{i});
        urt_adf_out.(Levelname).(clist{j}) = zeros(4,3);
        urt_ws_out.(Levelname).(clist{j}) = zeros(4,3);
        
        
        [junk posmin] = min(urt_adf.(Levelname_t).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(Levelname).(clist{j})(1,:) = [adfcrit1 urt_adf.(Levelname_t).(clist{j})(posmin,1) posmin];
        urt_ws_out.(Levelname).(clist{j})(1,:) = [wscrit1 urt_ws.(Levelname_t).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(Levelname_nt).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(Levelname).(clist{j})(2,:) = [adfcrit0 urt_adf.(Levelname_nt).(clist{j})(posmin,1) posmin];
        urt_ws_out.(Levelname).(clist{j})(2,:) = [wscrit0 urt_ws.(Levelname_nt).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(Idiffname).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(Levelname).(clist{j})(3,:) = [adfcrit0 urt_adf.(Idiffname).(clist{j})(posmin,1) posmin];
        urt_ws_out.(Levelname).(clist{j})(3,:) = [wscrit0 urt_ws.(Idiffname).(clist{j})(posmin,1) posmin];
        
        [junk posmin] = min(urt_adf.(IIdiffname).(clist{j})(:,info_urt)); %#ok
        urt_adf_out.(Levelname).(clist{j})(4,:) = [adfcrit0 urt_adf.(IIdiffname).(clist{j})(posmin,1) posmin]; 
        urt_ws_out.(Levelname).(clist{j})(4,:) = [wscrit0 urt_ws.(IIdiffname).(clist{j})(posmin,1) posmin];
               
    end
end




if not(gvnum==0)
    % global variables
    disp('- Testing global variables')
    for i=1:gvnum
        % initialize output structures, in which the following items will be stored: adf/ws stat,
        % aic, sbc, logl

        Levelname_t = sprintf('%s_t',gvnames{i});
        urt_adf.(Levelname_t) = zeros(maxlag_urt,4); % levels (with trend)
        urt_ws.(Levelname_t) = zeros(maxlag_urt,4);

        Levelname_nt = sprintf('%s_nt',gvnames{i});
        urt_adf.(Levelname_nt) = zeros(maxlag_urt,4); % levels (without trend)
        urt_ws.(Levelname_nt) = zeros(maxlag_urt,4);

        Idiffname = sprintf('D%s',gvnames{i});
        urt_adf.(Idiffname) = zeros(maxlag_urt,4); % Ist diff
        urt_ws.(Idiffname) = zeros(maxlag_urt,4);

        IIdiffname = sprintf('DD%s',gvnames{i});
        urt_adf.(IIdiffname) = zeros(maxlag_urt,4); % IInd diff
        urt_ws.(IIdiffname) = zeros(maxlag_urt,4);

        for p=1:maxlag_urt
            % testing levels
            % with trend
            urt_adf.(Levelname_t)(p,:) = adf(gv.(gvnames{i}),ct,p);
            urt_ws.(Levelname_t)(p,:) = ws(gv.(gvnames{i}),ct,p);
            % without trend
            urt_adf.(Levelname_nt)(p,:) = adf(gv.(gvnames{i}),c,p);
            urt_ws.(Levelname_nt)(p,:) = ws(gv.(gvnames{i}),c,p);
            % testing Ist diffs  (just intercept included)
            urt_adf.(Idiffname)(p,:) = adf(diff(gv.(gvnames{i}),1),c,p);
            urt_ws.(Idiffname)(p,:) = ws(diff(gv.(gvnames{i}),1),c,p);
            % testing IInd diffs (just intercept included)
            urt_adf.(IIdiffname)(p,:) = adf(diff(gv.(gvnames{i}),2),c,p);
            urt_ws.(IIdiffname)(p,:) = ws(diff(gv.(gvnames{i}),2),c,p);
        end

        % selecting ADF results depending on an information criterion
        urt_adf_out.(gvnames{i}) = zeros(4,3);
        urt_ws_out.(gvnames{i}) = zeros(4,3);


        [junk posmin] = min(urt_adf.(Levelname_t)(:,info_urt)); %#ok
        urt_adf_out.(gvnames{i})(1,:) = [adfcrit1 urt_adf.(Levelname_t)(posmin,1) posmin];
        urt_ws_out.(gvnames{i})(1,:) = [wscrit1 urt_ws.(Levelname_t)(posmin,1) posmin];

        [junk posmin] = min(urt_adf.(Levelname_nt)(:,info_urt)); %#ok
        urt_adf_out.(gvnames{i})(2,:) = [adfcrit0 urt_adf.(Levelname_nt)(posmin,1) posmin];
        urt_ws_out.(gvnames{i})(2,:) = [wscrit0 urt_ws.(Levelname_nt)(posmin,1) posmin];

        [junk posmin] = min(urt_adf.(Idiffname)(:,info_urt)); %#ok
        urt_adf_out.(gvnames{i})(3,:) = [adfcrit0 urt_adf.(Idiffname)(posmin,1) posmin];
        urt_ws_out.(gvnames{i})(3,:) = [wscrit0 urt_ws.(Idiffname)(posmin,1) posmin];

        [junk posmin] = min(urt_adf.(IIdiffname)(:,info_urt)); %#ok
        urt_adf_out.(gvnames{i})(4,:) = [adfcrit0 urt_adf.(IIdiffname)(posmin,1) posmin];
        urt_ws_out.(gvnames{i})(4,:) = [wscrit0 urt_ws.(IIdiffname)(posmin,1) posmin];

    end
end

