
function [wmat wmat_sy ylabelseq yearseq nyears cnames cnames_long cnum ...
    tvw_solution_modeflag] = build_wmat(intfname,fixedweights_flag,trade_tmp,...
    n_trdyears,trade_dates,firstobs,lastobs,min_trdyear,max_trdyear,cnum,...
    cnames,cnames_long,rnamesx,rnamesx_long,regionsx,regflag)

%**************************************************************************
% PURPOSE: retrieving cross-country (trade) flows and computing the weights 
% matrix           
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     

disp('- Retrieving flows data');
ccodes_t = xlsread(intfname,'MAIN','ccodes'); % import country codes
if isnumeric(ccodes_t(1))
    % let codes be string
    ccodes = {};
    for i=1:length(ccodes_t)
        ccodes = [ccodes; sprintf('%d',ccodes_t(i))]; %#ok
    end
end

allclist_t = trade_tmp(1,3:end);
% transform numeric to string
if isnumeric(allclist_t(1))
    allclist = {};
    for i=1:length(allclist_t)
        allclist = [allclist; sprintf('%d',allclist_t(i))]; %#ok
    end
end

dslct = []; % selection vector for country trade data
for j=1:length(allclist)
    for i=1:length(ccodes)
        if strcmp(allclist(j),ccodes(i))
            dslct = [dslct j]; %#ok
        end
    end
end

if fixedweights_flag == 1
    disp('- Fixed weight matrix will be computed, as selected');
    % empty matrices which refer to the time-varying weights case
    ylabelseq = [];
    yearseq = [];
    nyears = [];
    wmat_sy = [];
    tvw_solution_modeflag = [];

    first_trdyear = xlsread(intfname,'MAIN','first_trdyear');
    last_trdyear = xlsread(intfname,'MAIN','last_trdyear');

    for i=1:n_trdyears
        if trade_dates(i) == first_trdyear
            posfirst = i;
        end
        if trade_dates(i) == last_trdyear
            poslast = i;
        end
    end

    % import data on exports and imports of each country

    for n=1:cnum
        alltrd_tmp = xlsread('Flows\flows.xls',ccodes{n});

        alltrd = alltrd_tmp(2:end,3:end);
        % choose only data for countries existing in the current model
        trd = alltrd(:,dslct);
        % trim data according to period chosen
        trd_trm = trd(posfirst:poslast,:);
        if posfirst == poslast
            avg_trd = trd_trm;
        else
            % take averages over years
            avg_trd = mean(trd_trm);
        end
        % put data in structures
        tradedata.(cnames{n}) = avg_trd;
    end

    disp('- Building the weight matrix');
    trademat_tmp = [];
    tradematrix = zeros(cnum,cnum);
    for n=1:cnum
        trademat_tmp = [trademat_tmp; tradedata.(cnames{n})]; %#ok
    end
    % reorganize columns
    for j=1:cnum
        for i=1:cnum
            if isnan(trademat_tmp(i,j))
                tradematrix(:,i) = trademat_tmp(:,j);
            end
        end
    end

    % sobstitute NaNs with zeros
    tradematrix = nan2zero(tradematrix);
    % transpose the trade matrix
    tradematrix = tradematrix';

    if regflag == 1
        disp('- Updating weight matrix taking into account possible regions');
        % updating trade matrix in order to contain regions
        for r=1:length(rnamesx)
            [tradematrix cnames_long cnames cnum] = update_matrix(cnum,cnames,cnames_long,rnamesx(r),rnamesx_long(r),regionsx.(rnamesx{r}),tradematrix);
        end
    end

    % creating trade-based weights matrix
    wmat = weightmat(tradematrix);
    

elseif fixedweights_flag == 0 % time varying weights
    disp('- Time-varying weight matrices will be computed, as selected');
    
    
    trade_window = xlsread(intfname,'MAIN','trade_window');
    tvw_solution_modeflag = xlsread(intfname,'MAIN','tvw_solution_modeflag');
    
    firstobs_tmp = char(firstobs);
    firstyear_tmp = firstobs_tmp(1:4);
    firstyear = str2double(firstyear_tmp);

    lastobs_tmp = char(lastobs);
    lastyear_tmp = lastobs_tmp(1:4);
    lastyear = str2double(lastyear_tmp);

    nyears = lastyear-firstyear+1;
    yearseq = (firstyear:1:lastyear)';
    trdyearseq = (min_trdyear:1:max_trdyear)';

    if tvw_solution_modeflag == 1 % weights matrix for solution of the GVAR
        % is going to be constructed averaging over the
        % selected years window

        wmat_sy = [];

        % import data on exports and imports of each country
        for n=1:cnum

            alltrd_tmp = xlsread('Flows\flows.xls',ccodes{n});

            alltrd = alltrd_tmp(2:end,3:end);

            % choose only data for countries existing in the current model
            trd = alltrd(:,dslct);

%             % trim data according to period chosen
%             trd_trm = trd(posfirst:poslast,:);
            trd_trm = trd; % no need to trim trade data when doing t.v. weights. A.G. 27/10/2011
            

            ylabelseq = [];
            for t=1:nyears

                ylabel = sprintf('y%d',yearseq(t));
                ylabel_str = {ylabel};
                ylabelseq = [ylabelseq; ylabel_str]; %#ok

                if yearseq(t) <= min_trdyear + trade_window - 1
                    trd_tmp = trd_trm(1:trade_window,:);
                else
                    pointer = find(trdyearseq==yearseq(t));
                    trd_tmp = trd_trm(pointer-trade_window+1:pointer,:);
                end

                if trade_window > 1
                    % take averages over years if trade window is bigger than
                    % one year
                    avg_trd = mean(trd_tmp);
                else
                    avg_trd = trd_tmp;
                end

                % put data in structures
                tradedata.(cnames{n}).(ylabel) = avg_trd;
            end
        end

        disp('- Building the weight matrices for each sample year');

        wmat = [];      
        for t=1:nyears

            trademat_tmp = [];
            tradematrix = zeros(cnum,cnum);

            for n=1:cnum
                trademat_tmp = [trademat_tmp; tradedata.(cnames{n}).(ylabelseq{t})]; %#ok
            end

            % reorganize columns
            for j=1:cnum
                for i=1:cnum
                    if isnan(trademat_tmp(i,j))
                        tradematrix(:,i) = trademat_tmp(:,j);
                    end
                end
            end

            % sobstitute NaNs with zeros
            tradematrix = nan2zero(tradematrix);
            % transpose the trade matrix
            tradematrix = tradematrix';
            
            cnum_tmp = cnum;
            cnames_tmp = cnames;
            cnames_long_tmp = cnames_long;
            
            if regflag == 1
                if yearseq(t) == lastyear
                    disp('- Updating weight matrix taking into account possible regions');
                    % updating trade matrix in order to contain regions
                    for r=1:length(rnamesx)
                        [tradematrix cnames_long cnames cnum] = update_matrix(cnum,cnames,cnames_long,rnamesx(r),rnamesx_long(r),regionsx.(rnamesx{r}),tradematrix);
                    end
                else
                    for r=1:length(rnamesx)
                        [tradematrix cnames_long_tmp cnames_tmp cnum_tmp]  = update_matrix(cnum_tmp,cnames_tmp,cnames_long_tmp,rnamesx(r),rnamesx_long(r),regionsx.(rnamesx{r}),tradematrix);
                    end
                end
            end

            % creating trade-based weights matrix
            wmat(:,:,t) = weightmat(tradematrix); %#ok
        end

    else % two sets of weights matrices will be constructed: one in which
        % weights matrices are constructed using averages over the
        % years window, and another one of single-year flows weights
        % matrices

        % first let's create single-year flows weights matrices

        % backup of trade_window object
        trade_window_backup = trade_window;

        % set trade_window = 1
        trade_window = 1;


        for k=1:2
            % determine first position for trade-flows range
            if yearseq(1) <= trdyearseq(1) + trade_window -1
                posfirst = 1;
            else
                for i=1:length(trdyearseq)
                    if trdyearseq(i) == yearseq(1)
                        posfirst = i - trade_window + 1;
                        break
                    end
                end
            end

            % determine last position for trade-flows range
            if yearseq(end) <= trdyearseq(end)
                poslast = find(trdyearseq==yearseq(end));
            else
                error('Trade flows database is not up to date given the current estimation sample');
            end

            trdyearseq_trm = trdyearseq(posfirst:poslast);
            % import data on exports and imports of each country
            for n=1:cnum
                alltrd_tmp = xlsread('Flows\flows.xls',ccodes{n});
                alltrd = alltrd_tmp(2:end,3:end);
                % choose only data for countries existing in the current model
                trd = alltrd(:,dslct);
                % trim data according to period chosen
                trd_trm = trd(posfirst:poslast,:);
                                
                ylabelseq = [];
                for t=1:nyears
                    ylabel = sprintf('y%d',yearseq(t));
                    ylabel_str = {ylabel};
                    ylabelseq = [ylabelseq; ylabel_str]; %#ok

                    if yearseq(t) <= min_trdyear + trade_window - 1
                        trd_tmp = trd_trm(1:trade_window,:);
                    else      
                        pointer = find(trdyearseq_trm==yearseq(t));
                        trd_tmp = trd_trm(pointer-trade_window+1:pointer,:);
                    end

                    if trade_window > 1
                        % take averages over years if trade window is bigger than
                        % one year
                        avg_trd = mean(trd_tmp);
                    else
                        avg_trd = trd_tmp;
                    end
                    % put data in structures
                    tradedata.(cnames{n}).(ylabel) = avg_trd;
                end
            end


            if k==1
                disp('- Building the single-year weight matrices for each sample year');
            else
                disp('- Building the window-averaged weight matrices for each sample year');
            end

            wmat = [];
            for t=1:nyears

                trademat_tmp = [];
                tradematrix = zeros(cnum,cnum);

                for n=1:cnum
                    trademat_tmp = [trademat_tmp; tradedata.(cnames{n}).(ylabelseq{t})]; %#ok
                end

                % reorganize columns
                for j=1:cnum
                    for i=1:cnum
                        if isnan(trademat_tmp(i,j))
                            tradematrix(:,i) = trademat_tmp(:,j);
                        end
                    end
                end

                % sobstitute NaNs with zeros
                tradematrix = nan2zero(tradematrix);
                % transpose the trade matrix
                tradematrix = tradematrix';

                cnum_tmp = cnum;
                cnames_tmp = cnames;
                cnames_long_tmp = cnames_long;
                
                if regflag == 1
                    if k == 2 && yearseq(t) == lastyear
                        disp('- Updating weight matrix taking into account possible regions');
                        % updating trade matrix in order to contain regions
                        for r=1:length(rnamesx)
                            [tradematrix cnames_long cnames cnum] = update_matrix(cnum,cnames,cnames_long,rnamesx(r),rnamesx_long(r),regionsx.(rnamesx{r}),tradematrix);
                        end
                    else
                        for r=1:length(rnamesx)
                            [tradematrix cnames_long_tmp cnames_tmp cnum_tmp] = update_matrix(cnum_tmp,cnames_tmp,cnames_long_tmp,rnamesx(r),rnamesx_long(r),regionsx.(rnamesx{r}),tradematrix);
                        end
                    end
                end

                if k == 1 % creating single-year weights matrices
                    wmat_sy(:,:,t) = weightmat(tradematrix); %#ok
                else % creating years-window weights matrices
                    wmat(:,:,t) = weightmat(tradematrix); %#ok
                end
            end

            if k == 1 % restore original trade window
                trade_window = trade_window_backup;
            end
        end
    end
end