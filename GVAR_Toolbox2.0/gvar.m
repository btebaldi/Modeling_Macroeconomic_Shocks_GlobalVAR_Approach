
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all % clearing workspace
%**************************************************************************

disp('***************************************************************************');
disp('                                                                           ');
disp('                               GVAR Toolbox 2.0                            ');
disp('                                  August 2014                              '); 
disp('                                                                           ');
disp('                                                                           ');  
disp(' Alessandro Galesi, CEMFI Madrid & L.Vanessa Smith, University of York     ');
disp('***************************************************************************')


% Authors info:
%**************************************************************************
%
% Alessandro Galesig
% CEMFI, Madrid
% galesi@cemfi.edu.es
%
% L. Vanessa Smith
% University of York, UK
% vanessa.smith@york.ac.uk
%
%**************************************************************************


%% 0) INFO/TECH
% disp('********************************************************************');
format short g

addpath('Codes');


% program settings - do not modify
%***************************************
printout = 1;  % if set to 1, code prints output in output.xlsx
redoflag = 0; 
progdir = cd; 
%**************************************************************************

disp('>>> Type the interface filename (WITHOUT the .xls extension) that you would like to');
disp('    use and press enter. If you wish to use one of the demo interface files provided');
disp('    with the installation program, type gvarBriefDemo or gvarFullDemo, accordingly.');
disp('    Alternatively, type the name of your own interface file.');
disp('    Note that MatLab is case sensitive.');
disp(' ');
disp('======================================================================================');
disp('Note:');
disp('1. After typing the interface filename, the program will start running and it will make');
disp('   a number of pauses, unless the "Running the program with pauses option" is disabled.')
disp('   (This is the case for the gvarBriefDemo, for which no pauses will be performed). ');
disp('2. At each pause, you will be called upon to supply settings and/or check intermediate'); 
disp('   results in the MAIN worksheet of the interface file that will open automatically each');
disp('   time. Once this file has opened, always refer back to the MatLab command window for'); 
disp('   instructions and information.');
disp('3. Additional guidance is also available by clicking on many of the headings and field');
disp('   names within the MAIN worksheet.');
disp(' ');
disp('-----------------------------------------------------------------------------------');
disp('  After every pause, once the required settings have been supplied and/or the inter-');
disp('  mediate results have been checked, you must save and close the interface file. In');
disp('  fact, it is recommended that you close excel completely, each time.');
disp('-----------------------------------------------------------------------------------');
disp(' ');
disp(' USING THE FULL DEMO INTERFACE FILE');
disp('If you are using the full demo interface file most of the required settings and ');
disp('intermediate results are already provided. However, there are occassions where the');
disp('user is required to intervene. Once again, consult the MatLab command window for');
disp('instructions referring SPECIFICALLY to the use of the full demo interface file.');
disp('If no such instructions are given, just save and close the interface file whenever');
disp('it opens in order to proceed.');
disp('======================================================================================');
disp(' ');
intfname_tmp = input(' ','s');  
disp(' ');
disp('The program is now running (do not press any key)');



intfname = sprintf('%s.xls',intfname_tmp); 
if not(exist(intfname)) %#ok 
    % if there is no interface file with .xls extension, go for .xlsx 
    intfname = sprintf('%s.xlsx',intfname_tmp); % 
end
    
% special: check if pauses are allowed:
pauseflag = xlsread(intfname,'MAIN','pauseflag');  
% if set to 1, code features pauses
if isempty(pauseflag)
    pauseflag = 1;
end

%% 1) IMPORTING DATA
disp(' ');
disp('1) IMPORTING DATA');
disp('********************************************************************');


%% 1.1) Importing country and region names
%**************************************************************************

disp('1.1) Importing country and region names');

[junk cnames_long] = xlsread(intfname,'MAIN','cnames_long'); %#ok % countries names, format long
[junk cnames] = xlsread(intfname,'MAIN','cnames'); %#ok % countries names, format short 
cnum = length(cnames); % # of countries

% regions for estimation aggregation
[junk rnamesx_long] = xlsread(intfname,'MAIN','rnamesx_long'); %#ok % region names, format long 
[junk rnamesx] = xlsread(intfname,'MAIN','rnamesx'); %#ok % region names, format short 

if isempty(rnamesx_long) == 1
    regflag = 0; % no regional model aggregation
    regionsx = [];
    rnamesx_long = [];
    rnamesx = [];
else
    regflag =1;

    rnamesx_long = rnamesx_long(not(strcmp(rnamesx_long,'')));
    rnamesx = rnamesx(not(strcmp(rnamesx,'')));

    [junk regcountriesx] = xlsread(intfname,'MAIN','regcountriesx'); %#ok

    regcountriesx_temp = [];
    k=1;
    for i=1:length(regcountriesx)
        if strcmp(regcountriesx{i},'')
            regionsx.(rnamesx{k}) = regcountriesx_temp;
            k=k+1;
            regcountriesx_temp = [];
        else
            regcountriesx_temp = [regcountriesx_temp regcountriesx(i)]; %#ok
            if k==length(rnamesx)
                regionsx.(rnamesx{k}) = regcountriesx_temp;
            end
        end
    end

        clear junk regcountriesx_temp regcountriesx
end


%% 1.2) Importing weights for aggregation (used for regional analysis)
%**************************************************************************

aggrwgts_all = xlsread(intfname,'weights_aggr');
disp('1.2) Importing weights for aggregation (used for regional analysis)');
for n=1:cnum
    if not(isempty(aggrwgts_all))
        aggrwgts.(cnames{n}) = aggrwgts_all(n);
    else
        aggrwgts.(cnames{n}) = 0;
    end
end
clear aggrwgts_all

%% 1.3) Importing country data and GVAR settings
%**************************************************************************
disp('1.3) Importing the country specific datasets');

% import country-specific variable names
%***************************************
[junk vnames_long] = xlsread(intfname,'MAIN','vnames_long'); %#ok
[junk vnames] = xlsread(intfname,'MAIN','vnames'); %#ok
[junk vtypes_char] = xlsread(intfname,'MAIN','vtypes'); %#ok
vtypes = zeros(length(vtypes_char),1);
for i=1:length(vtypes)
    vtypec = vtypes_char(i);
    for j=1:10 % 10 is the max number of variable types
        vtname = sprintf('wmat%g',j);
        if strcmp(vtypec,vtname)
            vtypes(i) = j;
            break
        end
    end
end

    
vnames_long = vnames_long';
vnames = vnames';
vnum = length(vnames);

% Subsection dedicated to reading and printing first and last observational
% dates for the estimation window and for the cross-country flows
%**************************************************************************

[date_num date_chr] = xlsread(intfname,vnames{1},'A2:A65536');

if isempty(date_num)
    date = date_chr; % quarterly or monthly data
else
    date = date_num; % annual data
    % convert to character (fixed since version 1.2, A.G. 2012)
    date_chr = {};
    for i=1:length(date)
        datec = sprintf('%d',date(i));
        datec = {datec};
        date_chr = [date_chr; datec]; %#ok
    end
end

% Identifying the frequency of the data for use in forecasting

if not(iscell(date(1))); % annual data
    freq = 1;
else % data is quarterly or monthly
    date_cell = sprintf('%s',date{1});
    if  strcmp(date_cell(5),'Q'); % quarterly data
        freq = 4;
    elseif strcmp(date_cell(5),'M'); % monthly data
        freq = 12;
    end
end


% print initial obs and final obs in gvar.xls
min_date = date(1);
max_date = date(end);

% write the first and last obs on gvar.xls
xlswrite(intfname,[min_date max_date],'MAIN','Y15');

% import trade flows time range of the default database of the program

trade_tmp = xlsread('Flows\flows.xls');
trade_dates = trade_tmp(2:end,1);

min_trdyear = trade_dates(1);
max_trdyear = trade_dates(end);
n_trdyears = length(trade_dates);

% write the first and last trade year on gvar.xls
xlswrite(intfname,[min_trdyear max_trdyear],'MAIN','Y23');

% recover the number of weight types included in the interface file
[status,sheets] = xlsfinfo(intfname); % assuming gvar.xls is the name of the interface file
poswmat=find(strncmp(sheets, 'wmat', 4));
wnames=sheets(poswmat)';
xlswrite(intfname,wnames,'MAIN','A55');



if pauseflag == 1
% Ist pause: choose GVAR Settings
%**********************************************************************
winopen(intfname);
disp(' ');
msg = sprintf('>>> Pause and go to %s. Select the weights to be used for constructing',intfname);
disp(msg);  
disp('     the foreign variables in column N, and also define the GVAR settings for the:');
disp(' ');
disp('- Output folder name');
disp('- Frequency of pauses in the program');
disp('- Plotting of impulse responses');
disp('- Estimation window');
disp('- Type of weights');
disp('- Unit root tests');
disp('- Model selection');
disp('- Overidentifying restrictions test');
disp('- Weak exogeneity test');
disp('- Structural stability tests');
disp('- GVAR forecasts');
disp('- GVAR Trend/Cycle decomposition');
disp('- Dominant unit model');
disp(' ');
disp('  Once the above settings are defined (or if these are already defined),');
disp('  press enter.');
disp(' ');
disp('======================================================================================');
disp('When defining the settings for the program in the MAIN worksheet of the interface file:');
disp('1. White cells require you to fill in data, and blue cells contain information for you');
disp('   generated by the program.                                                       ');
disp('2. Dropdown menus are found in many of the fields. They are not apparent until you');
disp('   click within the field. Completing these fields is compulsory. For dropdown menus');
disp('   with the options 1 or 0: 1 means enable that function, and 0 disable that function.');
disp('3. If you choose to disable a function the program will ignore any other information');
disp('   relating to that function. This applies throughout the MAIN worksheet.');
disp('4. If you select one type of weights, it does not matter whether you leave the boxes');
disp('   related to alternative choices blank, as the program will ignore them by default.');
disp('   This applies throughout the MAIN worksheet of the interface file.                 ');
disp(' ');
disp(' USING THE FULL DEMO INTERFACE FILE');
disp('If you are using the full demo interface file and would like to replicate the results');
disp('in the Output Full Demo folder downloaded with the program, ensure that in column N');
disp('wmat1 is selected to construct the foreign variables of all domestic variables included');
disp('in the GVAR model.');
disp('======================================================================================');
disp(' ');
pause
msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
disp(msg);
endpause = input(' ','s'); 
disp('The program is now running (do not press any key)  ');
disp(' ');
end


%**************************************************************************
% the code now continues with the importing section
%**************************************************************************

%  Set output foldername
%**************************
[junk outdir_t] = xlsread(intfname,'MAIN','outdir'); %#ok
outdir_t = char(outdir_t);
outdir = sprintf('Output\\%s\\',outdir_t);
warning off all
mkdir(outdir);

% Set whether to plot graphs or not
%***********************************
graphsflag = xlsread(intfname,'MAIN','graphsflag'); % if =1, plots of irfs are performed

% Set modality of running the program
%*************************************
pauseflag = xlsread(intfname,'MAIN','pauseflag'); % if set to 1, code features pauses


%**************************************************************************
%                          Set estimation window
%**************************************************************************
[firstobs_num firstobs_chr] = xlsread(intfname,'MAIN','firstobs');
[lastobs_num lastobs_chr] = xlsread(intfname,'MAIN','lastobs');

if isempty(firstobs_num) % quarterly or monthly data
    annual=0;
    firstobs = firstobs_chr;
    lastobs = lastobs_chr;
    for i=1:length(date)
        if strcmp(date(i),firstobs)
            sposfirst = i;
        elseif strcmp(date(i),lastobs)
            sposlast = i;
        end
    end
    
    firstobs_tmp = char(firstobs);
    lastobs_tmp = char(lastobs);
    firstyear_estimation = str2num(firstobs_tmp(:,1:4)); %#ok
    lastyear_estimation = str2num(lastobs_tmp(:,1:4)); %#ok 
    
else % annual data
    annual=1;
    firstobs = firstobs_num;
    lastobs = lastobs_num;
    for i=1:length(date)
        if date(i)==firstobs
            sposfirst = i;
        elseif date(i)==lastobs
            sposlast = i;
        end
    end
    
    firstyear_estimation = firstobs;
    lastyear_estimation = lastobs;
end

nyears_estimation = lastyear_estimation-firstyear_estimation+1;
yearseq_estimation = firstyear_estimation:lastyear_estimation;

ylabelseq_estimation = [];
for t=1:nyears_estimation
    ylabel = sprintf('y%d',yearseq_estimation(t));
    ylabel_str = {ylabel};
    ylabelseq_estimation = [ylabelseq_estimation; ylabel_str]; %#ok
end

% Resize date according to sample estimation chosen
date = date(sposfirst:sposlast);

%***************************************************************************



% Unit root tests
urtflag = xlsread(intfname,'MAIN','urtflag');

% Structural stability tests
ss_flag = xlsread(intfname,'MAIN','ss_flag');

% Overidentifying restrictions tests
overid_flag = xlsread(intfname,'MAIN','overid_flag');

% Weak exogeneity test
we_flag = xlsread(intfname,'MAIN','we_flag');

% GVAR forecasts
% - unconditional forecasts
forc_flag = xlsread(intfname,'MAIN','forc_flag');
% - conditional forecasts
con_forc_flag = xlsread(intfname,'MAIN','con_forc_flag');


% Import country-specific variables
%*************************************
for i=1:vnum
    data_tmp = xlsread(intfname,vnames{i});
    if isequal(min_date,data_tmp(1,1)) && isequal(max_date,data_tmp(end,1))  % annual data case
        data_tmp = data_tmp(:,2:end);
    end
    %trimming
    data_trm = data_tmp(sposfirst:sposlast,:);
    alldata.(vnames{i}) = data_trm;
    
    if forc_flag == 1 && not(strcmp(lastobs,max_date))
        % saving also remaining actual data for eventually comparing it with GVAR
        % forecasts
        rdata_trm = data_tmp(sposlast+1:end,:);
        allrdata.(vnames{i}) = rdata_trm;
    end    
end

% Import global variables
%********************************
[junk gvnames_long] = xlsread(intfname,'MAIN','gvnames_long'); %#ok
[junk gvnames] = xlsread(intfname,'MAIN','gvnames'); %#ok

gvnames_long = gvnames_long';
gvnames = gvnames';
gvnum = length(gvnames);

if gvnum == 0
    gv = [];
    gvflag = [];
else
    for i=1:gvnum

        gv_tmp = xlsread(intfname,gvnames{i});
        if isequal(min_date,gv_tmp(1,1)) && isequal(max_date,gv_tmp(end,1)) % annual data case
            gv_tmp = gv_tmp(:,2:end);
        end        
        %trimming
        gv_trm = gv_tmp(sposfirst:sposlast,:);
        gv.(gvnames{i}) = gv_trm;
        
        if forc_flag == 1 && not(strcmp(lastobs,max_date))
            % saving also remaining actual data for eventually comparing it with GVAR
            % forecasts
            rgv_trm = gv_tmp(sposlast+1:end,:);
            rgv.(gvnames{i}) = rgv_trm;
        end           
    end
end


% Set the maximum lag order of the GVAR
%************************************************
maxlag_p = xlsread(intfname,'MAIN','maxlag_p');
maxlag_q = xlsread(intfname,'MAIN','maxlag_q');
maxlag = max(maxlag_p,maxlag_q);
% the maximum lag order of the GVAR is the maximum lag allowed for both 
% domestic and foreign variables across country-specific models


% Trend/Cycle decomposition of the GVAR model
TCflag = xlsread(intfname, 'MAIN','TCflag');% flag for doing Trend/Cycle decomposition

% Dominant unit model
dumodel_flag = xlsread(intfname,'MAIN','dumodel_flag');

% Feedback variables in the augmented regression
isfeedback = xlsread(intfname,'MAIN','addfeedbacks');



clear date_tmp junk data_tmp data_trm rdata_trm gv_tmp gv_trm rgv_trm


%% 1.4) Importing critical values for maximum eigenvalue and trace statistics
%**************************************************************************
disp('1.4) Importing critical values for maximum eigenvalue and trace statistics');
% note: critical values for cointegration tests with exogenous I(1)
% variables, at 95% significance and for case III and IV.
% See Mackinnon, Haug and Michelis (1999).

trace_crit95_c4_tmp = xlsread('Tech\coint_critvalues.xls','Trace_case4');
trace_crit95_c4 = trace_crit95_c4_tmp(2:end,2:end);
maxeig_crit95_c4_tmp = xlsread('Tech\coint_critvalues.xls','Maximum Eigenvalue_case4');
maxeig_crit95_c4 = maxeig_crit95_c4_tmp(2:end,2:end);
trace_crit95_c3_tmp = xlsread('Tech\coint_critvalues.xls','Trace_case3');
trace_crit95_c3 = trace_crit95_c3_tmp(2:end,2:end);
maxeig_crit95_c3_tmp = xlsread('Tech\coint_critvalues.xls','Maximum Eigenvalue_case3');
maxeig_crit95_c3 = maxeig_crit95_c3_tmp(2:end,2:end);
trace_crit95_c2_tmp = xlsread('Tech\coint_critvalues.xls','Trace_case2');
trace_crit95_c2 = trace_crit95_c2_tmp(2:end,2:end);
maxeig_crit95_c2_tmp = xlsread('Tech\coint_critvalues.xls','Maximum Eigenvalue_case2');
maxeig_crit95_c2 = maxeig_crit95_c2_tmp(2:end,2:end);

clear trace_crit95_c4_tmp maxeig_crit95_c4_tmp trace_crit95_c3_tmp maxeig_crit95_c3_tmp trace_crit95_c2_tmp maxeig_crit95_c2_tmp



%% 2) PREPARING DATA
disp(' ');
disp('2) PREPARING DATA');
disp('********************************************************************');



%% 2.1) Creating domestic variables for each country
%**************************************************************************
disp('2.1) Creating domestic variables for each country');

% creating structures for each variable
%****************************************
for i=1:vnum
    % determine the type of the variable
    dvtype.(vnames{i}) = vtypes(i);
    
    for n=1:cnum
        if alldata.(vnames{i})(:,n) == 123456789
            continue
        else
            dv.(vnames{i}).(cnames{n}) = alldata.(vnames{i})(:,n);
            if forc_flag == 1 && not(strcmp(lastobs,max_date))
                % saving also remaining actual data for eventually comparing it with GVAR
                % forecasts
                rdv.(vnames{i}).(cnames{n}) = allrdata.(vnames{i})(:,n);
            end        
        end
    end
end

nobs = length(dv.(vnames{1}).(cnames{1}));

if forc_flag == 1 && not(strcmp(lastobs,max_date))
    rnobs = length(rdv.(vnames{1}).(cnames{1}));
end
            

if regflag == 1
    %**************************************************************************
    % 2.1b) Creating region
    %**************************************************************************
    disp('2.1b) Creating regions ');

    for r=1:length(rnamesx)
        
        if annual == 0
            misal = not(strcmp(lastobs,max_date));
        else
            misal = not(isequal(lastobs,max_date));
        end
        
        if forc_flag == 1 && misal == 1
            % aggregating actual data not used for estimation, for comparing it with GVAR
            % forecasts
            rdv =  create_region(rnamesx(r), regionsx.(rnamesx{r}),rnobs, vnum, vnames, rdv, aggrwgts);
        end

        [dv aggrwgts] =  create_region(rnamesx(r), regionsx.(rnamesx{r}),nobs, vnum, vnames, dv, aggrwgts);
    end
    cnames_tmp = fieldnames(aggrwgts);
    cnum_tmp = length(cnames_tmp);
    cnames_long_tmp = [];
    for i = 1:cnum_tmp
        cname = cnames_tmp(i);
        rflag = 0;
        for k = 1:length(rnamesx)
            if strcmp(cname,rnamesx(k))
                cnames_long_tmp = [cnames_long_tmp; rnamesx_long(k)]; %#ok
                rflag = 1;
            end
        end
        if rflag == 0
            for j = 1:cnum
                if strcmp(cname,cnames(j))
                    cnames_long_tmp = [cnames_long_tmp; cnames_long(j)]; %#ok
                end
            end
        end
    end
else
    cnames_tmp = cnames;
    cnum_tmp = cnum;
    cnames_long_tmp = cnames_long;
end



% - domestic variables
%in the model's specification section
xlswrite(intfname,[vnames_long; vnames],'MAIN','AH3');
% in the simulation section
xlswrite(intfname,[vnames_long; vnames],'MAIN','ET3');

% - foreign-specific and feedback vars
xvec = []; fbvec = [];
for i=1:length(vnames)
    xvec{i} = [vnames{i} 's']; %#ok
    fbvec{i} = [vnames{i} '_tilde']; %#ok
end
xvec = xvec';
fbvec = fbvec';


%in the model's specification section
xlswrite(intfname,[vnames_long; xvec'],'MAIN','BC3');
%in the weak exogeneity test section
xlswrite(intfname,[vnames_long; xvec'],'MAIN','CT3');

% - global vars
if not(gvnum==0)
    % in the simulation section
    xlswrite(intfname,[gvnames_long; gvnames],'MAIN','FO3');
    %in the weak exogeneity test section
    xlswrite(intfname,[gvnames_long; gvnames],'MAIN','DO3');
    %in the model's specification section
    xlswrite(intfname,[gvnames_long; gvnames],'MAIN','BX3');
end

clear alldata allrdata


%% 2.2) Weight matrix
%**************************************************************************
disp('2.2) Weight matrix ');

% check if user wants to use different weight matrices for different variables
%**************************************************************************
multiplewmat = 0;
for i=1:length(vtypes)
    if max(vtypes) > 1
        multiplewmat = 1;
        break
    end
end
if multiplewmat == 1
    ntypes = max(vtypes);
    disp('- Multiple weight matrices are employed');
else
    ntypes = 1;
    disp('- Single weight matrix is employed');
end


% retrieve user's choice about weigths
%*************************************
[junk importbuild_wmat_flag]=xlsread(intfname,'MAIN','importbuild_wmat_flag'); %#ok
[junk weightstype_flag]=xlsread(intfname,'MAIN','weightstype_flag'); %#ok


% Weight matrix for type-1 variable (could be user-inputed or built by the program)
%**************************************************************************
if strcmp('program-built',importbuild_wmat_flag)
    disp('- Building the weight matrix using the flows.xls file');
    % the program builds the weight matrix using its own database
    % (it works just for variable type 1)
    
    if strcmp('fixed',weightstype_flag) == 1
        % fixed-weight matrix is chosen
        fixedweights_flag = 1;
    elseif strcmp('time-varying',weightstype_flag) == 1
        % time-varying weights are chosen
        fixedweights_flag = 0;
    else
        error('You have not selected an option for the type of weights');
    end
    [wmat wmat_sy ylabelseq yearseq nyears cnames cnames_long cnum tvw_solution_modeflag] = ...
        build_wmat(intfname,fixedweights_flag,trade_tmp,n_trdyears,trade_dates,firstobs,lastobs,min_trdyear,...
        max_trdyear,cnum,cnames,cnames_long,rnamesx,rnamesx_long,regionsx,regflag);   
    
    if printout == 1
        % print the weights matrix constructed from the toolbox
        if fixedweights_flag == 1
            disp('- Writing to output.xlsx: (fixed) weight matrix');
        elseif fixedweights_flag == 0
            disp('- Writing to output.xlsx: (time-varying) weight matrices');
        end
        print_weightmatrix(cnames_long,wmat,yearseq,weightstype_flag,outdir);
    end    
    
    % store in a structure the weights matrix
    if fixedweights_flag == 1
        % trick: generate a new matrix which vertically stacks T replicated weight matrices
        % where T is the number of years in the estimation sample
        wmat = repmat(wmat,nyears_estimation,1);
        wmat = reshape(wmat',cnum,cnum,nyears_estimation); % then reshape to get a 3d-matrix
        % then transpose weight matrices so that they add up to one by column
        for i=1:nyears_estimation
            wmat(:,:,i) = wmat(:,:,i)';
        end
    end
    wmatrices.wmat1 = wmat; % store the resulting sequence of matrices in a structure   
end

if (ntypes == 1 && strcmp('program-built',importbuild_wmat_flag))
    % no need to import additional weights matrices
else
    cnames = cnames_tmp;
    cnames_long = cnames_long_tmp;
    cnum = cnum_tmp;
    
    if strcmp('user-provided',importbuild_wmat_flag)
        start_type = 1; % type-1 variable has an user-inputed weight matrix
    elseif strcmp('program-built',importbuild_wmat_flag)
        start_type = 2;
    end
        
    disp('- Importing the weight matrix provided by the user');
    for k=start_type:ntypes
        % the program reads the weight matrix (or sequence of weight
        % matrices) from the interface file
        wmatsheetname = sprintf('wmat%d',k);
        [wdata labels] = xlsread(intfname,wmatsheetname);
        
        if length(wdata) == cnum
            % fixed-weight matrix
            %**********************
            wtype.(wmatsheetname) = 'f';
            
            wmat_tmp = wdata;
            wmat = [];
            for j=1:cols(wmat_tmp)
                if isnan(wmat_tmp(1,j))
                    continue
                else
                    wmat(:,j) = wmat_tmp(:,j); %#ok
                end
            end
            
            % generate a new matrix which vertically stacks T replicated weight matrices
            % where T is the number of years in the estimation sample
            wmat = repmat(wmat,nyears_estimation,1);     
        else
            % time-varying weight matrices
            %******************************
            wtype.(wmatsheetname) = 'tv';
            
            w_iniyear = wdata(1,1);
            w_lastyear = wdata(end,1);
            nyears = w_lastyear - w_iniyear + 1;
            yearseq = w_iniyear:w_lastyear;
            
            wmat = wdata(:,3:end);
        end
        wmat = reshape(wmat',cnum,cnum,nyears_estimation); % then reshape to get a 3d-matrix
        % then transpose weight matrices so that they add up to one by column
        for i=1:nyears_estimation
            wmat(:,:,i) = wmat(:,:,i)';
        end        
 
        wmatrices.(wmatsheetname) = wmat;
    end
end

% clear memory
clear trade_tmp


%% 3) COUNTRY MODELS

disp(' ');
disp('3) COUNTRY MODELS');
disp('********************************************************************');


%% 3.1 Model specification
%****************************************************************
disp('3.1) Model specification');

% create the array fvnames containing foreign variable names
%************************************************************
fvnames = [];
k=1;
for i=1:vnum
    fvnames{k} = vnames{i}; %#ok
    k=k+1;
end
fvnum = length(fvnames);

% for GIRFs & GFEVDs
xlswrite(intfname,[cnames_long cnames],'MAIN','EQ5');
% for model selection
xlswrite(intfname,[cnames_long cnames],'MAIN','AE5');
% check whether any specification is defined in the interface file
dvflag = xlsread(intfname,'MAIN','dvflag');


if isempty(dvflag) % the specification table is empty, force pauseflag=1
    pauseflag = 1;
end

if pauseflag == 1
    dvflag = ones(cnum,vnum);

    for j=1:vnum
        for i=1:cnum
            if isfield(dv.(vnames{j}),cnames(i))
            else
                dvflag(i,j) = NaN;
            end
        end
    end

    fvflag = ones(cnum,fvnum);
    if not(gvnum==0)
        gvflag = ones(cnum,gvnum);
    else
        gvflag = [];
    end
    
    
    % for model selection
    if not(gvnum==0)
        xlswrite(intfname,gvflag,'MAIN','BX5');
    end
    xlswrite(intfname,fvflag,'MAIN','BC5');
    xlswrite(intfname,dvflag,'MAIN','AH5');
    
    
    % pause and let user modify model specification
    %************************************************
    winopen(intfname);
    disp(' ');
    msg = sprintf('>>> Pause and go to %s: Define the specification of the individual models,',intfname);
    disp(msg);
    disp('   and then press enter.');
    disp(' ');
    disp('=========================================================================================');
    disp('1. By default the program generates the value of 1 for all variables (domestic, foreign');
    disp('   and global) in all country models, that have been defined at the outset of the MAIN');
    disp('   worksheet. A value of 1 indicates inclusion of that variable in the model. To exclude');
    disp('   a variable from the model, set its value to 0. To include a global variable as domestic');
    disp('   set its value to 2.                     ');
    disp('2. If a dominant unit model is NOT enabled, every global variable has to enter as endogenous');
    disp('   in ONE country only, and can enter as weakly exogenous in any or all remaining countries.');
    disp('3. If a dominant unit model is enabled, any global variable that will subsequently enter');  
    disp('   this model should NOT be included as a domestic variable in the GVAR model i.e. no ');
    disp('   value of 2 should appear anywhere in the column for this global variable. It can enter');
    disp('   as weakly exogenous in any or all countries. Those global variables that will not enter');
    disp('   the dominant unit model should be treated as under point 2 above.');
    disp('4. DO NOT fill any empty cells. An empty cell means that the data is not available for that');
    disp('   particular variable in that country.                                          ');
    disp('5. For every run of the program the specification of the individual models needs to be');
    disp('   redefined. To retain the specification in a subsequent run, disable the "Run the');
    disp('   program with pauses" function at the initial settings stage, and select no from the');
    disp('   dropdown list in the lag order selection field of the model selection display.');
    disp(' ');
    disp(' USING THE FULL DEMO INTERFACE FILE');
    disp('If you are using the full demo interface file, to define the specification of the');
    disp('country-specific models which follows largely DdPS(2007) except for the oil price variable:');
    disp('1. Exclude the domestic variable, ep, from the US model (set it to 0).');
    disp('2. Exclude the foreign variable, eps, from all country models (set it to 0) except the US.');
    disp('3. Exclude the foreign variables, eqs, rs and ls, from the US model (set them to 0).');
    disp('4. Include all global variables as weakly exogenous in all country models. DO NOT include');
    disp('   any of the global variables as domestic in any country model.'); 
    disp('=========================================================================================');
    disp(' ');
    pause
    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
    disp(msg);
    endpause = input(' ','s');
    disp('The program is now running (do not press any key)  ');
    disp(' ');
end


% extract info about country specification
%******************************************
dvflag = xlsread(intfname,'MAIN','dvflag');
fvflag = xlsread(intfname,'MAIN','fvflag');
if not(gvnum==0)
    gvflag = xlsread(intfname,'MAIN','gvflag');
end
    

% Checking consistency of GVAR specification
%************************************************************************
doloop = 0; % loop for checking which global variables are included as endogenous in the GVAR
weloop = 0; % loop for checking that the number of weakly exogenous variables does
% not exceed 8 for each individual model (when doing cointegration test)

while doloop < 1 || weloop < 1
    
    % Check that the number of weak exogenous variables does not exceed 8.
    %***********************************************************************
    dvflag = xlsread(intfname,'MAIN','dvflag');
    fvflag = xlsread(intfname,'MAIN','fvflag');
    if not(gvnum==0)
        gvflag = xlsread(intfname,'MAIN','gvflag');
    end
    
    
    rank_tmp =xlsread(intfname,'MAIN','rank'); % to check if cells of rank are empty
    if isempty(rank_tmp) || pauseflag == 1
        % the code is going to run the cointegration test
        
        if not(gvnum==0)
            sumwevars = sum(fvflag,2) + sum(gvflag==1,2);
        else
            sumwevars = sum(fvflag,2);
        end
        
        if any(sumwevars>8)
            disp('');
            disp('Warning: In each individual model the sum of foreign and global variables cannot exceed eight.');
            disp('');
            
            % pause and let the user change the specification of the model
            %************************************************************
            winopen(intfname);
            disp(' ');
            msg = sprintf('>>>> Pause and go to %s: Change the specification of the models so that each model includes ',intfname);
            disp(msg);
            disp('        at most eight foreign and global variables.');
            pause
            disp(' ');
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s');
            disp('The program is now running (do not press any key)  ');
            disp(' ');
        else
            weloop = 1;
        end
    else
        weloop = 1;
    end
    
    
    % Check which variables enter in the GVAR as global
    %***********************************************************
    
    gxidx = zeros(gvnum,1);
    
    % global (endogenous) variables
    gnv = [];
    gnvnames = [];
    gnvnames_long = [];
    
    % global (exogenous) variables
    gxv = [];
    gxvnames = [];
    gxvnames_long = [];
    
    isdumodel = 0;
    if not(gvnum==0)
        for j=1:gvnum
            flag = 0;
            for i=1:rows(gvflag)
                if gvflag(i,j) == 2
                    flag = 1; % global variable is endogenous
                    gnv.(gvnames{j}) = gv.(gvnames{j});
                    gnvnames = [gnvnames gvnames(j)]; %#ok
                    gnvnames_long = [gnvnames_long gvnames_long(j)]; %#ok
                end
            end
            if flag == 0 % global variable is exogenous --> a dominant unit model exists
                isdumodel = 1;
                
                gxidx(j) = 1;
                
                gxv = [gxv gv.(gvnames{j})]; %#ok
                gxvnames = [gxvnames gvnames(j)]; %#ok
                gxvnames_long = [gxvnames_long gvnames_long(j)]; %#ok
            end
        end
    end
    
    gxvnum = cols(gxv);
    gnvnum = gvnum - gxvnum;
    
    
    duerror = 0;
    if dumodel_flag == 1 && isdumodel == 0
        % the user has specified that a dominant unit model exists but then all
        % global variables are endogenous
        duerror = 1;
    elseif dumodel_flag == 0 && isdumodel == 1
        % the user has specified that a dominant unit model does not exist but
        % then one or more global variables are exogenous
        duerror = 2;
    end
    
    
    if duerror > 0;
        winopen(intfname);
        if duerror == 1
            disp('>>> Warning: Make sure that at least one global variable is not included as domestic');
            disp('    (endogenous) in any of the individual country models (i.e. ensure that no value');
            disp('     of "2" appears in the column of at least one global variable).');
        elseif duerror == 2
            disp('>>> Warning: Make sure you have defined each global variable to be domestic in ONE individual model.');
        end
        disp(' ');
        pause
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
        
        gvflag = xlsread(intfname,'MAIN','gvflag');
    else
        % all is ok
        doloop = 1;
    end
end



% delete domestic variables which are de-selected in the model specification section:
%**************************************************************************
for j=1:vnum
    for i=1:cnum
        if dvflag(i,j) == 0
            dv.(vnames{j}) = rmfield(dv.(vnames{j}),cnames{i});
            if forc_flag == 1 && not(strcmp(lastobs,max_date))
               rdv.(vnames{j}) = rmfield(rdv.(vnames{j}),cnames{i});
            end
        end
    end
end



%% 3.2) Creating foreign-specific variables
%**************************************************************************
disp('3.2) Creating foreign-specific variables ');


% generate structures of domestic, foreign variables and weights (where
% weights could be either fixed or time-varying)
if freq == 1 % annual data
    ydate = date;
else
    char_date = char(date);
    ydate = str2num(char_date(:,1:4)); %#ok
end


for t=1:nyears_estimation
    
    % generate the structure of domestic variables
    obs_idx = [];
    for j=1:nobs
        if  ydate(j) == yearseq_estimation(t)
            obs_idx = [obs_idx; j]; %#ok
        end
    end
    for i=1:vnum
        clist = fieldnames(dv.(vnames{i}));
        for n=1:length(clist)
            dv_t.(ylabelseq_estimation{t}).(vnames{i}).(clist{n}) = dv.(vnames{i}).(clist{n})(obs_idx);
        end
    end
    
    % generate the structure of weight matrices
    for k=1:ntypes
        wmattype = sprintf('wmat%d',k);    
        wm_t.(ylabelseq_estimation{t}).(wmattype) = wmatrices.(wmattype)(:,:,t);
    end
    
    % generate the structure of foreign variables
    fv_t.(ylabelseq_estimation{t}) = create_foreignvariables(cnum, cnames, vnum, vnames, dv_t.(ylabelseq_estimation{t}),dvtype, wm_t.(ylabelseq_estimation{t}));
end

% now assemble the series of foreign-specific variables
for t=1:nyears_estimation
    if t==1
        fv = struct(fv_t.(ylabelseq_estimation{t}));
    else
        % reconstruct the series:
        for i=1:vnum
            clist = fieldnames(fv_t.(ylabelseq_estimation{t}).(vnames{i}));
            for n=1:length(clist)
                fv.(vnames{i}).(clist{n}) = [fv.(vnames{i}).(clist{n}); fv_t.(ylabelseq_estimation{t}).(vnames{i}).(clist{n})];
            end
        end
    end
end

% ### cleaning variables from memory ###
clear fv_t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if printout == 1
    % Descriptive stats of domestic, foreign-specific and global variables
    disp('- Adding to output.xlsx: Descriptive statistics of domestic, foreign-specific and global variables');
    print_dstats(vnum,vnames,vnames_long,cnum,cnames,cnames_long,dv,...
        fv,gv,gvnum,gvnames,gvnames_long,outdir);
end



if urtflag == 1
    %% 3.3) Unit root tests
    %****************************************************************
    disp('3.3) Unit root tests');
    
    maxlag_urt = xlsread(intfname,'MAIN','urt_maxlag');
    [junk info_urt_choice] = xlsread(intfname,'MAIN','urt_criterion'); %#ok
    

    [urt_adf_out urt_ws_out] = unitroot_tests(vnum,fvnum,gvnum,vnames,gvnames,dv,fv,gv,maxlag_urt,info_urt_choice);


    if printout == 1
        disp('- Adding to output.xlsx: Unit root tests');
        print_unitroottests(vnum,vnames,cnum,cnames,cnames_long,gvnum,gvnames,urt_adf_out,urt_ws_out,outdir);
    end

    clear urt_adf_out urt_adf_stats urt_ws_out urt_ws_stats
end


%% 3.4) Building country models
%**************************************************************************
disp('3.4) Building country models');


[endog endoglist junk exog exoglist exognx exognxlist] = create_countrymodels(cnames,cnum, ...
    vnames,vnames_long,fvnames,vnum,gvnames,gvnames_long,gvnum,gnvnum,gnvnames,dv,fv,gv,dvflag,fvflag,gvflag); %#ok


clear endoglist_long

if annual == 0
    misal = not(strcmp(lastobs,max_date));
else
    misal = not(isequal(lastobs,max_date));
end

if forc_flag == 1 && misal == 1
    % saving remaining actual data for comparing it with GVAR
    % forecasts

    for n=1:cnum
        rendog.(cnames{n}) = [];   % creating the block of endogenous variables
        for i=1:vnum
            if dvflag(n,i) == 1
                rendog.(cnames{n}) = [rendog.(cnames{n}) rdv.(vnames{i}).(cnames{n})];
            end
        end

        if not(gvnum==0)
            % adding the global variables
            for i=1:gvnum
                if gvflag(n,i) == 2 % global variable is endogenous
                    rendog.(cnames{n}) = [rendog.(cnames{n}) rgv.(gvnames{i})];
                end
            end
        end
    end
end


if printout == 1
    
    if annual == 1 % if using annual, use dates in string format
        date = date_chr;
    end
    
    disp('- Writing to countrydata.xlsx: Data for domestic and foreign variables of each country');
    print_countrydata(date,cnum,cnames,cnames_long,...
        endoglist,endog,exoglist,exog,gvnames,outdir);
    
    if dumodel_flag == 1
    % print also the dominant unit model data
    disp('- Adding to countrydata.xlsx: Data for dominant unit model');
    vlab = date;
    hlab = {'date'};
    tab = [];
    for i=1:length(gxvnames)
        hlab = [hlab gxvnames(i)]; %#ok
        tab = [tab gxv(:,i)]; %#ok
    end
    toprint = [hlab; vlab num2cell(tab)];
    xlswrite([outdir 'countrydata.xlsx'],toprint,'DOMINANT UNIT');  
    end
end



%% 3.5) Determining the lag length of each country model
%**************************************************************************
disp('3.5) Determining the lag orders of each country model');


% retrieve lag order criterion
%*****************************
[junk criterion] = xlsread(intfname,'MAIN','criterion'); %#ok

if strcmp(criterion,'aic') == 1
    lagselect = 1;
    info = 2;
elseif strcmp(criterion,'sbc') == 1
    lagselect = 1;
    info = 3;
elseif strcmp(criterion,'no') == 1
    lagselect = 0;
end
% Information criterion used for model lag selection:
% info = 2 for Akaike Information criterion;
% info = 3 for Schwartz Bayesian Information Criterion.

% retrieve maximum lag order for serial correlation test
%*******************************************************
psc = xlsread(intfname,'MAIN','psc');

if lagselect == 0  % user provides arbitrary lag orders
    
    varxlag_tmp= xlsread(intfname,'MAIN','varxlag');
    
else % lag order selection is performed

    for n=1:cnum
        aic_sbc.(cnames{n}) = [];
        F_serialcorr.(cnames{n}) = [];
        for p=1:maxlag_p
            for q=1:maxlag_q
                [logl_aic_sbc  psc_degfrsc_Fcrit_Fsc] = select_varxlag(maxlag,psc,endog.(cnames{n}),p,exog.(cnames{n}),q);
                aic_sbc.(cnames{n}) = [aic_sbc.(cnames{n}); logl_aic_sbc p q];
                F_serialcorr.(cnames{n}) = [F_serialcorr.(cnames{n}); psc_degfrsc_Fcrit_Fsc];
            end
        end
    end
    for n=1:cnum
        varxlag.(cnames{n}) = [];
        for i=1:rows(aic_sbc.(cnames{1}))
            if aic_sbc.(cnames{n})(i,info) == max(aic_sbc.(cnames{n})(:,info))
                varxlag.(cnames{n}) = aic_sbc.(cnames{n})(i,4:5);
            end
        end
    end
    
    varxlag_tmp = [];
    for n=1:cnum
        varxlag_tmp = [varxlag_tmp; varxlag.(cnames{n})]; %#ok
    end
end

clear psc_degfrsc_Fcrit_Fsc

% print to the interface file the resulting lag orders
%**************************************************************************
xlswrite(intfname,varxlag_tmp,'MAIN','CJ5');

if pauseflag == 1
    % pause and let user eventually modify lag orders
    %************************************************
    winopen(intfname);
    disp(' ');
    msg = sprintf('>>> Pause and go to %s: Check the lag orders found (or inputted if no selection criterion',intfname);
    disp(msg);   
    disp('    was previously used), then press enter.');
    pause
    disp(' ');
    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
    disp(msg);
    endpause = input(' ','s');
    disp('The program is now running (do not press any key)  ');
    disp(' ');

    % retrieve lag orders from gvar.xls
    varxlag_tmp= xlsread(intfname,'MAIN','varxlag');
end



for n=1:cnum
    varxlag.(cnames{n}) = varxlag_tmp(n,:);
end

clear varxlag_tmp

if printout == 1  % if lag order selection was performed
    if lagselect == 1
        disp('- Adding to output.xlsx: VARX* order selection results and residual serial correlation test results');

        % print results in 'output.xlsx'
        print_VARX_SC(cnum,cnames,cnames_long,aic_sbc,F_serialcorr,vnames,gvnames,vnum,gvnum,dvflag,gvflag,outdir);
        
        clear F_serialcorr aic_sbc        
    end

    disp('- Adding to output.xlsx: VARX* lag orders');
    % print VARX* lag orders in 'output.xlsx'
    print_VARXord(cnames,cnames_long,cnum,varxlag,outdir);
end

clear fsc


%% 3.6) Determining the number of cointegrating relations for each country
%**************************************************************************
disp('3.6) Determining the number of cointegrating relations for each country');

estcase_tmp =xlsread(intfname,'MAIN','estcase');

if isempty(estcase_tmp) || pauseflag==1
    estcase_tmp = [];
    for n=1:cnum
        estcase_tmp = [estcase_tmp; 4]; %#ok % write on gvar.xls the default option of estimation case IV
    end
    xlswrite(intfname,estcase_tmp,'MAIN','CN5');
end

if pauseflag == 1
    % pause and let user choose estimation case
    %************************************************
    winopen(intfname);
    disp(' ');
    msg = sprintf('>>> Pause and go to %s: Select between case II (restricted intercept),',intfname);
    disp(msg);
    disp('     case III (unrestricted intercept) and case IV (unrestricted intercept, restricted trend)');
    disp('     for VECMX* estimation, then press enter.');
    pause
    disp(' ');
    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
    disp(msg);
    endpause = input(' ','s');
    disp('The program is now running (do not press any key)  ');
    disp(' ');
end

% retrieve VECMX* estimation cases for each country-model
%******************************************************
estcase_tmp =xlsread(intfname,'MAIN','estcase');

for n=1:cnum
    estcase.(cnames{n}) = estcase_tmp(n);
end


% compute cointegration statistics
%*********************************

rank_tmp =xlsread(intfname,'MAIN','rank');

if isempty(rank_tmp) || pauseflag == 1
    for n=1:cnum
        [junk trace_out maxeig_out] = cointegration_test(maxlag,estcase.(cnames{n}),endog.(cnames{n}),...
            varxlag.(cnames{n})(1),exog.(cnames{n}), varxlag.(cnames{n})(2)); %#ok

        trace.(cnames{n}) = trace_out;
        maxeig.(cnames{n}) = maxeig_out;
    end


    % determining the number of cointegrating relations
    %**************************************************
    [rank trace_critvals] = get_rank(cnum,cnames,endog,exog,trace,trace_crit95_c2,trace_crit95_c3,...
        trace_crit95_c4,maxeig,maxeig_crit95_c2,maxeig_crit95_c3,maxeig_crit95_c4,estcase);

    clear lambda 

    % print ranks found in gvar.xls
    %**************************************************************************
    rank_tmp = [];
    for n=1:cnum
        rank_tmp = [rank_tmp; rank.(cnames{n})]; %#ok
        if rank.(cnames{n}) == length(endoglist.(cnames{n}))
            msg = sprintf('>>> Warning: Full rank found for model of %s',cnames_long{n});
            disp(msg);
        end
    end
    xlswrite(intfname,rank_tmp,'MAIN','CQ5');


    if printout == 1
        disp('- Adding to output.xlsx: Cointegration test statistics');
        % Printing cointegration statistics to output.xlsx
        print_cointmaxtraceVARX(cnum,cnames,cnames_long,endoglist,exoglist,maxeig,trace,trace_critvals,outdir);

        clear maxeig trace trace_critvals
    end
end

if pauseflag == 1
    % pause and let user eventually modify ranks
    %**********************************************************************
    %****
    winopen(intfname);
    disp(' ');
    msg = sprintf('>>> Pause and go to %s: Check the ranks found, then press enter.',intfname);
    disp(msg);
    disp(' ');
    pause
    disp(' ');
    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
    disp(msg);
    endpause = input(' ','s');
    disp('The program is now running (do not press any key)  ');
    disp(' ');
end

% retrieve ranks from gvar.xls
%**************************************************************************
rank_tmp = xlsread(intfname,'MAIN','rank');

for n=1:cnum
    rank.(cnames{n}) = rank_tmp(n);
end

if printout == 1
    disp('- Adding to output.xlsx: Cointegrating ranks');
    % printing number of cointegrating relations
    %*******************************************
    print_ncntVARX(cnum,cnames,cnames_long,rank,outdir);
end

clear rank_tmp


if overid_flag == 1
    %% Impose overidentifying restrictions on the cointegrating relations 
    %**************************************************************************
    disp(' ');
    disp('- Imposing overidentifying restrictions on the cointegrating relations');
 
    [beta_r nunrestrpar] = overid_restr(cnum,cnames,cnames_long,endoglist,...
         exoglist,gvnum,gvnames,estcase,outdir,rank);


    if isempty(beta_r)
        overid_flag = 0;
    end
else
    beta_r = [];
end


%% 3.7) Estimate individual VECMX*
%**************************************************************************
disp('3.7) Estimating individual VECMX* models ');

% Estimate VECMX* for each country (exactly identified estimation)
%****************************************************************
for n = 1: cnum
    [beta.(cnames{n}) alpha.(cnames{n}) Psi.(cnames{n}) ...
        epsilon.(cnames{n}) Omega.(cnames{n}) ecm.(cnames{n}) ...
        std.(cnames{n}) logl.(cnames{n}) aic.(cnames{n}) ...
        sbc.(cnames{n}) r2.(cnames{n}) rbar2.(cnames{n}) ...
        aics.(cnames{n}) sbcs.(cnames{n}) hcwstd.(cnames{n}) ...
        nwcstd.(cnames{n}) psc_degfrsc_Fcrit_Fsc.(cnames{n})] = mlcoint(maxlag,rank.(cnames{n}), ...
        estcase.(cnames{n}),psc,endog.(cnames{n}),varxlag.(cnames{n})(1), ...
        exog.(cnames{n}),varxlag.(cnames{n})(2));
end



if overid_flag == 1
    % Perform overidentifying restriction estimation
    %******************************************
        
    for n=1:cnum
        if isfield(beta_r,cnames(n))
            
            [alpha.(cnames{n}) Psi.(cnames{n}) epsilon.(cnames{n}) ...
                Omega.(cnames{n}) ecm.(cnames{n}) std.(cnames{n}) ...
                logl_r.(cnames{n}) aic.(cnames{n}) sbc.(cnames{n}) r2.(cnames{n}) ...
                rbar2.(cnames{n}) aics.(cnames{n}) sbcs.(cnames{n}) ...
                hcwstd.(cnames{n}) nwcstd.(cnames{n}) ...
                psc_degfrsc_Fcrit_Fsc.(cnames{n})] = mlcoint_r(beta_r.(cnames{n}),...
                maxlag,estcase.(cnames{n}),endog.(cnames{n}),varxlag.(cnames{n})(1),...
                exog.(cnames{n}),varxlag.(cnames{n})(2),psc);
            
            beta.(cnames{n}) = beta_r.(cnames{n});
            
            overid_LR.(cnames{n}) = -2*(logl_r.(cnames{n})-logl.(cnames{n}));
            overid_dgf.(cnames{n}) = rows(beta.(cnames{n}))*cols(beta.(cnames{n})) ...
                                    - (cols(beta.(cnames{n})))^2 - nunrestrpar.(cnames{n});
        else
            overid_LR.(cnames{n}) = NaN;
            overid_dgf.(cnames{n}) = NaN;
        end
    end
end



% Normalize beta vectors (of exactly identified estimation)
%*********************************************************
for n=1:cnum
    if not(isfield(beta_r,cnames(n)))
        rk = rank.(cnames{n});
        Pi_tmp = alpha.(cnames{n}) *beta.(cnames{n})';
        
        if estcase.(cnames{n}) == 4 || estcase.(cnames{n}) == 2
            Pi = Pi_tmp(:,2:end); % remove trend, if estcase = 4; intercept, if estcase = 2
            alpha_tmp = Pi(:,1:rk);
            alpha_norm.(cnames{n}) = alpha_tmp;
            
            invb = eye(rows(beta.(cnames{n})(2:rk+1,1:rk)))/beta.(cnames{n})(2:rk+1,1:rk);
        elseif estcase.(cnames{n}) == 3
            Pi = Pi_tmp; % no trend
            alpha_tmp = Pi(:,1:rk);
            alpha_norm.(cnames{n}) = alpha_tmp;
            
            invb = eye(rows(beta.(cnames{n})(1:rk,1:rk)))/beta.(cnames{n})(1:rk,1:rk);
        end
        beta_norm.(cnames{n}) = beta.(cnames{n})*invb;
    else
        alpha_norm.(cnames{n}) = alpha.(cnames{n});
        beta_norm.(cnames{n}) = beta_r.(cnames{n});
    end
end


if printout == 1
    disp('- Adding to output.xlsx: VECMX* estimates');
    % printing VECMX* estimates
    %************************
    print_ECMS_VARX(cnum,cnames,cnames_long,gvnames,endoglist,exoglist,varxlag,Psi,alpha,beta,estcase,outdir);

    disp('- Adding to output.xlsx: VECMX* statistics');
    % printing VECMX* statistics
    %*************************
    tab = ['VECMX* Estimation: Statistics' num2cell(NaN(1,3))];
    tab = [tab; num2cell(NaN(1,4))];
    tab = [tab; {'Country' 'logLik' 'Akaike' 'Schwartz Bayesian'}];

    for n=1:cnum
        tab = [tab; cnames_long(n) num2cell(logl.(cnames{n})) ...
            num2cell(aic.(cnames{n})) num2cell(sbc.(cnames{n}))]; %#ok
    end
    
    xlswrite([outdir 'output.xlsx'],tab,'ECMS_stats');

    clear aic sbc

    disp('- Adding to output.xlsx: VECMX* cointegrating vectors');
    % printing VECMX* cointegrating vectors
    %*************************************
    print_ECMS_CVs(cnum,cnames,endoglist,exoglist,cnames_long,rank,alpha_norm,beta_norm,estcase,gvnames,outdir);

    disp('- Adding to output.xlsx: VECMX* single-equation statistics');

    % printing VECMX* single-equation statistics
    %*************************************
    tab = ['VECMX* Single-Equation Statistics' num2cell(NaN(1,5))];
    tab = [tab; num2cell(NaN(1,6))];
    tab = [tab; {'Country' 'Variable' 'R_square' 'Rbar_square' 'AIC' 'SBC'}];
    for n=1:cnum
        vlist = endoglist.(cnames{n});
        for i=1:length(vlist)
            stat = [r2.(cnames{n})(i) rbar2.(cnames{n})(i) aics.(cnames{n})(i) sbcs.(cnames{n})(i)]; 
            tab = [tab; cnames_long(n) vlist(i) num2cell(stat)]; %#ok
        end
    end
    
    xlswrite([outdir 'output.xlsx'],tab,'ECMS_seq_stats');
    
    clear title hlab vlab tab aics rbar2 r2 sbcs

    disp('- Adding to output.xlsx: Descriptive statistics and normality test of VECMX* residuals');
    % printing descriptive statistics (& normality tests) of VECMX* residuals
    %**********************************************************************
    print_ECMS_resstats(cnum,cnames,cnames_long,endoglist,epsilon,outdir);
    
    disp('- Adding to output.xlsx: F-statistics for the test of serial correlation of the VECMX* residuals');
    % printing serial correlation test statistics of VECMX* residuals
    %**********************************************************************
    print_ECMS_restestsc(cnum,cnames,cnames_long,psc_degfrsc_Fcrit_Fsc,...
        vnames,gvnames,vnum,gvnum,dvflag,gvflag,outdir);
  
end


clear Pi Pi_tmp alpha_tmp invb Omega alpha_norm psc_degfrsc_Fcrit_Fsc


if we_flag == 1
    %% 3.8) Testing Weak Exogeneity
    %**************************************************************************
    disp('3.8) Testing weak exogeneity  ');
    
    % retrieve lag order criterion
    %*****************************
    [junk criterion_we] = xlsread(intfname,'MAIN','criterion_we'); %#ok
    
    if strcmp(criterion_we,'aic') == 1
        lagselect_we = 1;
        info_we = 2;
    elseif strcmp(criterion_we,'sbc') == 1
        lagselect_we = 1;
        info_we = 3;
    elseif strcmp(criterion_we,'no') == 1
        lagselect_we = 0;
    end
    % Information criterion used for weak exogeneity test
    % info = 2 for Akaike Information criterion;
    % info = 3 for Schwartz Bayesian Information Criterion.
    
    
    % read cells of lag orders
    % ************************************************************
    varxlag_we_tmp= xlsread(intfname,'MAIN','varxlag_we');    
    
    if isempty(varxlag_we_tmp)
        pauseflag = 1;
    end
    
    if  pauseflag==1 
        % print default settings for regressors (the ones used in estimation)
        spec_we_tmp = [];
        gvflag_we_tmp = [];
        
        if not(gvnum == 0)
            for j=1:cols(gvflag)
                for i=1:rows(gvflag)
                    if gvflag(i,j) == 1
                        gvflag_we_tmp(i,j) = 1; %#ok
                    elseif gvflag(i,j) == 0
                        gvflag_we_tmp(i,j) = 0 ; %#ok
                    elseif gvflag(i,j) == 2
                        gvflag_we_tmp(i,j) = 0; %#ok
                    end
                end
            end
        end
        
        fvflag_we_tmp = fvflag;
        
        if not(gvnum == 0)
            xlswrite(intfname,gvflag_we_tmp,'MAIN','DO5');
        end
        xlswrite(intfname,fvflag_we_tmp,'MAIN','CT5');
        
        
        % pause and let user eventually modify model specification
        %**********************************************************************
        %****
        winopen(intfname);
        disp(' ');
        msg = sprintf('>>> Pause and go to %s: Check the regression specifications for the weak exogeneity test,',intfname);
        disp(msg);
        disp('    and then press enter.');
        disp(' ');
        disp('=========================================================================================');
        disp('Note that for these regression specifications:');
        disp('1. By default, the program will set the same lag orders and same specification for the');
        disp('   domestic, foreign and global variables as used at the estimation stage of the country-');
        disp('   specific models.                                                                      ');
        disp('2. If you wish, you can change the lag orders and specification for the foreign and global');
        disp('   variables. The specification for the domestic variables remains set by default.');
        disp('3. Any global variable previously included as domestic in the estimation stage of a country-');
        disp('   specific model will automatically be excluded (set to 0) from the specification of that');
        disp('   model, when testing for weak exogeneity. This should NOT be changed, as for each country');
        disp('   the domestic variables are included in the the weak exogeneity test regressions by default.');
        disp('4. For every run of the program changes in the regression specification need to be re-entered.');
        disp('   To retain the specification in a subsequent run, disable the "Run the program with pauses"');
        disp('   function at the initial settings stage.');
        disp(' ');
        disp(' USING THE FULL DEMO INTERFACE FILE');
        disp('If you are using the full demo interface file to define the regression specification for the');
        disp('weak exogeneity test, as in DdPS(2007):                    ');
        disp('1. Include the foreign variable, eps, in all country models by setting it to 1.')
        disp('=========================================================================================')
        disp(' ');
        pause
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
    end
        
    % retrieve specification of foreign specific vars for weak exog tests
    %******************************************************************
    fvflag_we = xlsread(intfname,'MAIN','fvflag_we');
    
    if not(gvnum == 0)
        gvflag_we_tmp = xlsread(intfname,'MAIN','gvflag_we');
        gvflag_we = [];
        for j=1:cols(gvflag_we_tmp)
            if isnan(gvflag_we_tmp(1,j))
                continue
            else
                gvflag_we(:,j) = gvflag_we_tmp(:,j); %#ok
            end
        end
    end
    
    for n=1:cnum
        exog_we.(cnames{n}) = [];  % creating the block of (weakly) exogenous variables
        exoglist_we.(cnames{n}) = {}; % creating the list of exogenous var names
        
        for i=1:fvnum
            if fvflag_we(n,i) == 1
                exog_we.(cnames{n}) = [exog_we.(cnames{n}) fv.(fvnames{i}).(cnames{n})];
                exoglist_we.(cnames{n}) = [exoglist_we.(cnames{n}) fvnames{i}];
            end
        end
        
        if not(gvnum == 0)
            for j=1:gvnum
                if gvflag_we(n,j) == 1
                    exog_we.(cnames{n}) = [exog_we.(cnames{n}) gv.(gvnames{j})];
                    exoglist_we.(cnames{n}) = [exoglist_we.(cnames{n}) gvnames{j}];
                end
            end
        end
    end        
    

    % Lag orders of weak exogeneity test
    % **********************************
    if  lagselect_we == 1 % lag order selection is performed

        % retrieve maximum lag order for serial correlation test
        %*******************************************************
        psc_we = xlsread(intfname,'MAIN','psc_we');
    
        % retrieve maximum lag orders of WE regression equations
        %*******************************************************
        px_we = xlsread(intfname,'MAIN','px_we'); 
        qx_we = xlsread(intfname,'MAIN','qx_we');
        
        for n=1:cnum
            aic_sbc_we.(cnames{n}) = [];
            F_serialcorr_we.(cnames{n}) = [];
            for p=1:px_we
                for q=1:qx_we
                   [logl_aic_sbc_we  psc_degfrsc_Fcrit_Fsc_we] = select_lags_we(psc_we,...
                    endog.(cnames{n}),p,exog.(cnames{n}),q,exog_we.(cnames{n}),ecm.(cnames{n})); 
                                       
                   aic_sbc_we.(cnames{n}) = [aic_sbc_we.(cnames{n}); logl_aic_sbc_we p q];
                   F_serialcorr_we.(cnames{n}) = [F_serialcorr_we.(cnames{n}); psc_degfrsc_Fcrit_Fsc_we];
                end
            end
        end
        for n=1:cnum
            varxlag_we.(cnames{n}) = [];
            for i=1:rows(aic_sbc_we.(cnames{1}))
                if aic_sbc_we.(cnames{n})(i,info_we) == max(aic_sbc_we.(cnames{n})(:,info_we))
                    varxlag_we.(cnames{n}) = aic_sbc_we.(cnames{n})(i,4:5);
                end
            end
        end
        
        varxlag_we_tmp = [];
        for n=1:cnum
            varxlag_we_tmp = [varxlag_we_tmp; varxlag_we.(cnames{n})]; %#ok
        end
        
    % print to the interface file the resulting lag orders
    %**************************************************************************
    xlswrite(intfname,varxlag_we_tmp,'MAIN','DZ5');   
    
    end


    if pauseflag == 1
        
        if lagselect_we == 0
            % print default settings for lag orders (the ones used in estimation)
            varxlag_we_tmp = [];
            for n=1:cnum
                varxlag_we_tmp = [varxlag_we_tmp; varxlag.(cnames{n})]; %#ok
            end
            xlswrite(intfname,varxlag_we_tmp,'MAIN','DZ5');
        end        
        
        % pause and let user eventually modify lag orders
        %************************************************
        winopen(intfname);
        disp(' ');
        msg = sprintf('>>> Pause and go to %s: Check the lag orders found (or imposed if no selection criterion',intfname);
        disp(msg);
        disp('    was previously used) for the weak exogeneity regressions, then press enter.');
        pause
        disp(' ');
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
    end
    
    % retrieve lag orders from gvar.xls
    varxlag_we_tmp= xlsread(intfname,'MAIN','varxlag_we');  
    for n=1:cnum
        varxlag_we.(cnames{n}) = varxlag_we_tmp(n,:);
    end
    clear varxlag_we_tmp
    
    
    
    % do weak exogeneity test
    %**************************************************************************
    wetest_stat = test_weakexogeneity(cnames,cnum,endog,exog,varxlag_we,exog_we,ecm,rank);
    
    
    if printout == 1
        if lagselect_we == 1  % if lag order selection was performed 
            disp('- Adding to output.xlsx: Lag order selection results for the weak exogeneity regressions');
            disp('  and residual serial correlation test results');
            print_orderwereg_sc(cnum,cnames,cnames_long,aic_sbc_we,F_serialcorr_we,vnames,gvnames,vnum,gvnum,fvflag,gvflag,outdir)

            clear F_serialcorr_we aic_sbc_we
        end        
        
        disp('- Adding to output.xlsx: Weak exogeneity test results');
        print_exogeneity_test(cnum,cnames,cnames_long,varxlag_we,fvnames,gvnames,wetest_stat,exoglist,outdir);     
    end
    
    
    clear exog_we exoglist_we varxlag_we fvflag_we gvflag_we gvflag_we_tmp wetest_stat
end



%% 3.9) Contemporaneous Effects of Foreign Variables on their Domestic Counterparts
%**************************************************************************
disp('3.9) Contemporaneous effects of foreign variables on their domestic counterparts');

% retrieve contemporaneous coefficients
%**************************************
[ctpcoeffs ctpstd ctpstd_hcw ctpstd_nwc ctptvals ...
    ctptvals_hcw ctptvals_nwc ctplabel]=contmpcoeff(cnum,...
    vnames,endoglist,exoglist,cnames,Psi,std,hcwstd,nwcstd,estcase);

if printout == 1
    disp('- Adding to output.xlsx: Contemporaneous coefficients');
    % print contemporaneous coefficients in output.xlsx
    %*************************************************
    print_ContmpCoeff(cnum,cnames,cnames_long,ctplabel,ctpcoeffs,ctpstd,...
        ctpstd_hcw,ctpstd_nwc,ctptvals,ctptvals_hcw,ctptvals_nwc,vnames,outdir);
end

% clearing
clear ctplabel ctpcoeffs ctpstd ctpstd_hcw ctpstd_nwc ctptvals ctptvals_hcw ctptvals_nwc ...
    std hcwstd nwcstd

%% 3.10) Average pairwise cross-section correlations: Variables and residuals
%**************************************************************************
disp('3.10) Average pairwise cross-section correlations: Variables and residuals');


% compute average pairwise cross section correlations
[avgcorr avgcorr_d avgcorr_VARXres] = avgcorrs(cnum,cnames,maxlag,dv,endoglist,vnames,epsilon);

if printout == 1
    disp('- Adding to output.xlsx: Average pairwise cross-section correlations');
    % print avgcorrs
    print_avgcorr_vall(cnum,cnames_long,vnames,avgcorr,avgcorr_d,avgcorr_VARXres,outdir);

    clear avgcorr avgcorr_d avgcorr_VARXres
end



if ss_flag == 1
    %% 3.11) Structural stability tests
    %**************************************************************************
    disp('3.11) Structural stability tests');


    ccut_tmp = xlsread(intfname,'MAIN','chow_cut'); % Percentage of trimming for the structural stability chow test for the country-specific VARX* models @
    ccut = ccut_tmp/100;

    [kpsup kpmnsq ny rny qlr mw apw rqlr rmw rapw maxobs] = structural_stability_tests(cnum,cnames,endoglist,endog,exog,varxlag,estcase,maxlag,ecm,ccut);

    if printout == 1
        disp('- Adding to output.xlsx: Structural stability test statistics');
        print_sstests(kpsup,kpmnsq,ny,rny,qlr,mw,apw,rqlr,rmw,rapw,maxobs,...
            cnum,cnames,cnames_long,date,maxlag,outdir);
    end

    clear ccut_tmp ecm maxobs

end




if gxvnum > 0 
    %%                   Model for the dominant unit
    %********************************************************************** 
    if ss_flag == 1
    disp('3.12) Specification and estimation of the dominant unit model');
    else
    disp('3.11) Specification and estimation of the dominant unit model');    
    end
      
    % print info on DOMINANT UNIT sheet
    %*********************************************************************
    if gxvnum == 1
        xlswrite(intfname,{'Univariate'},'DOMINANT UNIT','G5');
    else
        xlswrite(intfname,{'Multivariate'},'DOMINANT UNIT','G5');
    end
    
    
    % print names of global variables
    xlswrite(intfname,[gxvnames_long' gxvnames'],'DOMINANT UNIT','M5');
    
    if isfeedback == 1
        % print names of feedback variables
        xlswrite(intfname,[vnames_long; fbvec'],'DOMINANT UNIT','S3');
        xlswrite(intfname,fbvec,'DOMINANT UNIT','AG5');
        %print names of countries associated with the weights to construct
        %feedback variables
        xlswrite(intfname,cnames_long,'DOMINANT UNIT','AK5');
        
        if pauseflag == 1
            % print selection matrix for feedback variables
            xlswrite(intfname,ones(gxvnum,vnum),'DOMINANT UNIT','S5');
        end
    end
    

    if pauseflag == 1
    % pause and let user modify model specification
    %************************************************
    winopen(intfname);
    disp(' ');
    msg = sprintf('>>> Pause and go to %s: Define the settings and the specification of the',intfname);
    disp(msg);
    disp('    dominant unit model. This includes information required in the columns adjacent to');
    disp('    the settings (if applicable), in particular:');
    disp('   - the maximum/actual lag order for the global variables');
    disp('   - the specification of the feedback variables');
    disp('   - the maximum/actual lag order for the feedback variables');
    disp('   - the weights to construct the feedback variables');
    disp('   - specifying which weights from the drop down lists (column AH)');
    disp('     to be used for computing the feedback variables that will be');
    disp('     included in the dominant unit model,')
    disp('   then press enter.');
    disp(' ');
        disp(' USING THE FULL DEMO INTERFACE FILE');
        disp('If you are using the full demo interface file and would like to replicate the results');
        disp('in the Output Full Demo folder:                                                      ');
        disp(' 1. Ensure that the value of 2 is defined as the maximum lag order for both the VEC model');
        disp('    and the augmented dominant unit model (for all equations).');
        disp(' 2. Include only the y_tilde and Dp_tilde feedback variables in the augmented regression');
        disp('    by changing the 1s to 0s for the rest of the feedback variables.');
        disp(' 3. Ensure that in column AH wmat1 is selected to construct the above two feedback variables.');
        disp('======================================================================================== ');
    disp(' ');
    pause
    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
    disp(msg);
    endpause = input(' ','s');
    disp('The program is now running (do not press any key)  ');
    disp(' ');
    end
    
    % if isfeedback = 1, check whether the user has inputed the feedback variables
    %**********************************************************************
    if isfeedback == 1
        doloop = 0;
        while doloop < 1
            
            feedbacks_flagmat = xlsread(intfname,'DOMINANT UNIT','feedbacks');
            fbvars_byeq = sum(feedbacks_flagmat,2);
            
            duerror = 0;
            sum_fbvars_byeq = sum(fbvars_byeq);
            if sum_fbvars_byeq == 0
                duerror = 1;
            end
            isfeedback_byeq = zeros(length(fbvars_byeq),1);
            for i=1:length(fbvars_byeq)
                if fbvars_byeq(i)>0
                    isfeedback_byeq(i) = 1;
                end
            end
            

            if duerror > 0;
                winopen(intfname);
                if gxvnum == 1
                    disp('>>> Warning: Make sure you have defined at least one feedback variable in the augmented');
                    disp('    univariate model.');
                elseif gxvnum > 1
                    disp('>>> Warning: Make sure you have defined at least one feedback variable in at least one equation');
                    disp('    of the augmented VEC model.');
                end
                disp(' ');
                pause
                msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                disp(msg);
                endpause = input(' ','s');
                disp('The program is now running (do not press any key)  ');
                disp(' ');
                
                gvflag = xlsread(intfname,'MAIN','gvflag');
            else
                % all is ok
                doloop = 1;
            end
        end
    end

    % Retrieve info on specification of dominant unit model
    %**********************************************************************
    
    if gxvnum > 1 % Multivariate model
        
        % retrieve the maximum lag order of the dominant unit model
        maxlag_du = xlsread(intfname,'DOMINANT UNIT','maxlag_du');
        
        % retrieve lag order criterion
        [junk criterion_du] = xlsread(intfname,'DOMINANT UNIT','criterion_du'); %#ok
        if strcmp(criterion_du,'aic') == 1
            lagselect_du1 = 1;
            info = 2;
        elseif strcmp(criterion_du,'sbc') == 1
            lagselect_du1 = 1;
            info = 3;
        elseif strcmp(criterion_du,'no') == 1
            lagselect_du1 = 0;
        end
        % Information criterion used for model lag selection:
        % info = 2 for Akaike Information criterion;
        % info = 3 for Schwartz Bayesian Information Criterion.
        
        % retrieve maximum lag order for serial correlation test
        psc_du = xlsread(intfname,'DOMINANT UNIT','psc_du');
        
        if lagselect_du1 == 0  % user provides arbitrary lag orders
            
            varxlag_du = maxlag_du; % the maximum lag order is the actual lag order
            aic_sbc_du = [];
            F_serialcorr_du = [];
            
        else % lag order selection is performed
            if pauseflag == 1 % do lag order selection
                aic_sbc_du = [];
                F_serialcorr_du = [];
                for p=1:maxlag_du
                    
                    exog_du = []; q = 0; % VEC model
                    [logl_aic_sbc  psc_degfrsc_Fcrit_Fsc] = select_varxlag(maxlag_du,psc_du,gxv,p,exog_du,q);
                    
                    aic_sbc_du = [aic_sbc_du; logl_aic_sbc p q]; %#ok
                    F_serialcorr_du = [F_serialcorr_du; psc_degfrsc_Fcrit_Fsc]; %#ok
                end
                
                varxlag_du = [];
                for i=1:rows(aic_sbc_du)
                    if aic_sbc_du(i,info) == max(aic_sbc_du(:,info))
                        varxlag_du = aic_sbc_du(i,4:5);
                    end
                end
            else
                % retrieve the lag order previously found by lag order
                % selection criteria
                varxlag_du = xlsread(intfname,'DOMINANT UNIT','varxlag_du');
            end
        end
        
        
        
        if pauseflag == 1
            if not(lagselect_du1 == 0) % lag order selection has been chosen
                % print lag order found in gvar.xls
                xlswrite(intfname,varxlag_du(1),'DOMINANT UNIT','J14');
            end
        end
        
        if not(lagselect_du1 == 0) % lag order selection has been chosen
            if pauseflag == 1
                % pause and let user eventually modify the lag order for the dominant unit model
                %************************************************
                winopen(intfname);
                disp(' ');
                msg = sprintf('>>> Pause and go to %s: Check the lag order found for the dominant unit VAR model',intfname);
                disp(msg);
                disp('     (or inputted if no selection criterion was previously used), then press enter.');
                pause
                disp(' ');
                msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                disp(msg);
                endpause = input(' ','s');
                disp('The program is now running (do not press any key)  ');
                disp(' ');
            end
            % retrieve lag order from gvar.xls from varxlag_du, both with
            % and without pauses
            varxlag_du= xlsread(intfname,'DOMINANT UNIT','varxlag_du');
        end
        
        
        % retrieve the case of treatment of deterministics in VECM
        [junk estcase_du_char] = xlsread(intfname,'DOMINANT UNIT','deterministics_vecm_du'); %#ok
        if strcmp(estcase_du_char,'Case 4');
            estcase_du = 4;
        elseif strcmp(estcase_du_char,'Case 3');
            estcase_du = 3;
        elseif strcmp(estcase_du_char,'Case 2');
            estcase_du = 2;
        end
        
        
        % Estimate the VECM model (Ist stage)
        %************************************
        
        % Compute cointegration statistics
        stat_du = 0; % stat_du: 0 for trace, 1 for max eigenvalue; use trace statistic as default
        
        rank_du =xlsread(intfname,'DOMINANT UNIT','rank_du');
        
        if isempty(rank_du) || pauseflag == 1
            
            [lambda_out trace_out maxeig_out] = cointegration_test(maxlag_du,estcase_du,gxv,varxlag_du(1));
            
            lambda_du = lambda_out;
            trace_du = trace_out;
            maxeig_du = maxeig_out;
            
            % Determining the number of cointegrating relations
            rank_du =0;
            
            if estcase_du == 4
                
                critvals_tmp = trace_crit95_c4(1:gxvnum,1); % note: here use the critical values of column k=0
                trace_critvals_du = flipud(critvals_tmp);
                
                critvals_tmp = maxeig_crit95_c4(1:gxvnum,1); % note: here use the critical values of column k=0
                maxeig_critvals_du = flipud(critvals_tmp);
                
            elseif estcase_du == 3
                
                critvals_tmp = trace_crit95_c3(1:gxvnum,1); % note: here use the critical values of column k=0
                trace_critvals_du = flipud(critvals_tmp);
                
                critvals_tmp = maxeig_crit95_c3(1:gxvnum,1); % note: here use the critical values of column k=0
                maxeig_critvals_du = flipud(critvals_tmp);
                
            elseif estcase_du == 2
                
                critvals_tmp = trace_crit95_c2(1:gxvnum,1); % note: here use the critical values of column k=0
                trace_critvals_du = flipud(critvals_tmp);
                
                critvals_tmp = maxeig_crit95_c2(1:gxvnum,1); % note: here use the critical values of column k=0
                maxeig_critvals_du = flipud(critvals_tmp);
            end
            
            if stat_du == 0  % Trace statistic chosen
                for i=1:gxvnum
                    if trace_du(i) > trace_critvals_du(i)
                        rank_du=rank_du+1;
                    else
                        break
                    end
                end
            elseif stat_du == 1  % Maximum eigenvalue statistic chosen
                for i=1:gxvnum
                    if maxeig_du(i) > maxeig_critvals_du(i)
                        rank_du=rank_du+1;
                    else
                        break
                    end
                end
            end
            
            if rank_du == gxvnum
                msg = sprintf('>>> Warning: Full rank found in model for the dominant unit');
                disp(msg);
            end
            
            % print rank found in gvar.xls
            xlswrite(intfname,rank_du,'DOMINANT UNIT','G21');
            
            if pauseflag == 1
                if printout == 1
                    disp('- Adding to output.xlsx: VAR order selection and cointegration results for the dominant unit VEC model');
                    print_DUmodel_VECM_stats(gxvnum,gxvnames,aic_sbc_du,F_serialcorr_du,varxlag_du(1),rank_du,maxeig_du,trace_du,...
                        trace_critvals_du,outdir);
                end
            end
            
            
            if pauseflag == 1
                % pause and let user eventually modify the rank of the dominant
                % unit model
                %******************************************************************
                winopen(intfname);
                disp(' ');
                msg = sprintf('>>> Pause and go to %s: Check the rank found for the dominant unit VEC model, then press enter.',intfname);
                disp(msg);
                disp(' ');
                pause
                disp(' ');
                msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                disp(msg);
                endpause = input(' ','s');
                disp('The program is now running (do not press any key)  ');
                disp(' ');
            end
            
            % retrieve rank from gvar.xls
            rank_du = xlsread(intfname,'DOMINANT UNIT','rank_du');
        end
        
        % VECM estimation
        % skip if there are no coefficients to estimate
        if (rank_du == 0) && (estcase_du == 2) && (varxlag_du(1)  == 1)
            % there are no coefficients to estimate in the VEC model
            
            ecm_du = [];
            alpha_du = [];
            beta_du = [];
            beta_norm_du = [];
            beta_norm.du_model = [];
            
        else
            [beta_du alpha_du Psi_du epsilon Omega ecm_du std logl aic sbc r2 rbar2 ...
                aics sbcs hcwstd nwcstd psc_degfrsc_Fcrit_Fsc_du] = mlcoint(maxlag_du,...
                rank_du,estcase_du,psc_du,gxv,varxlag_du(1));
            
            
            % Normalize the beta vector
            if estcase_du == 4 || estcase_du == 2 % skip trend if case=4, intercept if case=2
                invb = eye(rows(beta_du(2:rank_du+1,1:rank_du)))/beta_du(2:rank_du+1,1:rank_du);
            elseif estcase_du == 3 % no trend
                invb = eye(rows(beta_du(1:rank_du,1:rank_du)))/beta_du(1:rank_du,1:rank_du);
            end
            beta_norm_du = beta_du*invb;
            
            % store beta_du and beta_norm_du in the beta and beta_norm
            % structures (which contain the beta vectors of country models)
            beta.du_model = beta_du;
            beta_norm.du_model = beta_norm_du;
            % also store the treatment of the deterministic component in the
            % VECM of the dominant unit model, as it will be useful later in
            % case that the beta vector has to be reordered
            estcase.du_model = estcase_du;
            % also store the cointegrating rank of the dominant unit model
            rank.du_model = rank_du;
            
            if printout == 1
                disp('- Adding to output.xlsx: Estimates of the dominant unit model (VEC model)');
                print_DUmodel_VECM_est(gxvnames,varxlag_du,estcase_du,beta_du,alpha_du,Psi_du,psc_degfrsc_Fcrit_Fsc_du,outdir);
            end
        end
    else
        maxlag_du = [];
        varxlag_du = [];
        rank_du = [];
        ecm_du = [];
        alpha_du = [];
        beta_du = [];
        beta_norm_du = [];
        beta_norm.du_model = [];
        estcase_du = [];
    end
    
    
    % Retrieve info about the augmented regression
    %***********************************************
    
    % read max (or actual) lag orders of global variables
    max_ptildel = xlsread(intfname,'DOMINANT UNIT','max_ptildel');
    
    % retrieve lag order criterion
    [junk criterion_du_2] = xlsread(intfname,'DOMINANT UNIT','criterion_du_2'); %#ok
    if strcmp(criterion_du_2,'aic') == 1
        lagselect_du2 = 1;
        info2 = 2;
    elseif strcmp(criterion_du_2,'sbc') == 1
        lagselect_du2 = 1;
        info2 = 3;
    elseif strcmp(criterion_du_2,'no') == 1
        lagselect_du2 = 0;
        info2 = [];
    end
    % Information criterion used for model lag selection:
    % info = 2 for Akaike Information criterion;
    % info = 3 for Schwartz Bayesian Information Criterion.
    
    % read whether there are feedbacks
    feedbacks_flagmat = xlsread(intfname,'DOMINANT UNIT','feedbacks');    
    
    
    if isfeedback == 1
        k = 1;
        while k < 2
            % read choices and data for weights used to compute the feedback variables
            % and check that every feedback variable has its own
            % associated set of weights
            errflag = 0;
            
            [junk fbvarstypes] = xlsread(intfname,'DOMINANT UNIT','fwvars'); %#ok
            feedback_dataweights = xlsread(intfname,'DOMINANT UNIT','feedback_dataweights');
            
            if not(cols(fbvarstypes)>1)
                errflag = 1; % the user has not defined any weight
            else
                
                fbtypes_char = fbvarstypes(:,2);
                fbnames = fbvarstypes(:,1);
                fbnum = length(fbvarstypes(:,1));
                
                fbtypes = zeros(length(fbtypes_char),1);
                for i=1:length(fbtypes)
                    fbtypec = fbtypes_char(i);
                    for j=1:6 % 6 is the max number of variable types for the feedbacks
                        vtname = sprintf('wmat%g',j);
                        if strcmp(fbtypec,vtname)
                            fbtypes(i) = j;
                            break
                        end
                    end
                end
                
                % check that every feedback variable has its own associated set of weights
                if rows(feedbacks_flagmat) == 1 % univariate model
                    isfbvaridx = feedbacks_flagmat;
                elseif rows(feedbacks_flagmat) > 1 % multivariate model
                    isfbvaridx = max(feedbacks_flagmat);
                end
                
                
                for i=1:length(isfbvaridx)
                    if isfbvaridx(i) > 0 && fbtypes(i) == 0
                        % a feedback variable has not its own associated set of
                        % weights
                        errflag = 1;
                    end
                end
            end
            
            if errflag == 0
                k = 2; % all is ok
            else
                disp('>>> Warning: Missing associated weight matrix for some feedback variables.')
                disp(' ');
                msg = sprintf('>>> Pause and go to %s: Ensure that each feedback variable of those included',intfname);
                disp(msg);
                disp('  in the dominant unit model has an associated weight matrix, then press enter.');                              
                disp(' ');
                winopen(intfname);
                pause
                msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                disp(msg);
                disp('The program is now running (do not press any key)  ');
                disp(' ');
            end
        end
        
        
        for n=1:cnum
            feedwgts.(cnames{n}) = feedback_dataweights(n,:);
        end
        
        % calculating weights
        fweights = country_weights(feedwgts,cnames,vnum,vnames,gnvnum,gnvnames,endoglist,fbtypes);       
                                            
        % build feedback variables
        for g = 1:gxvnum
            
            exog_du.(gxvnames{g}) = [];
            exognum_du.(gxvnames{g}) = sum(feedbacks_flagmat(g,:));
            
            % identify which feedback variables will enter in equation of
            % global variable g
            fvdunames.(gxvnames{g}) = {};
            for k=1:length(vnames)
                if feedbacks_flagmat(g,k) == 1
                    fvdunames.(gxvnames{g}) = [fvdunames.(gxvnames{g}) vnames(k)];
                end
            end
            
            % create the feedback variables
            fvlist = fvdunames.(gxvnames{g});
            for j = 1 : length(fvlist)
                clist = fieldnames(dv.(fvlist{j}));
                
                sumx = 0;
                
                for i = 1:length(clist)
                    sumx = sumx + fweights.(fvlist{j}).(clist{i})*dv.(fvlist{j}).(clist{i});
                end
                
                fvdu.(fvlist{j}) = sumx;
                exog_du.(gxvnames{g}) = [exog_du.(gxvnames{g}) fvdu.(fvlist{j})];
            end
            
            exognames_du.(gxvnames{g}) = fvdunames.(gxvnames{g});
        end
        
        % read max (or actual) lag orders of feedback variables
        max_qtildel_tmp = xlsread(intfname,'DOMINANT UNIT','max_qtildel');
        max_qtildel = max_qtildel_tmp(:,2);
        if length(max_qtildel) < gxvnum
            max_qtildel = [max_qtildel; NaN(gxvnum-length(max_qtildel),1)];
        end
        
        for i=1:length(isfeedback_byeq)
            if isfeedback_byeq(i) == 0
                max_qtildel(i) = 1; % conventional order (feedback variables won't be included in estimation)
            end
        end
                
    else
        for g = 1:gxvnum
            exog_du.(gxvnames{g}) = [];
            exognum_du.(gxvnames{g}) = [];
            exognames_du.(gxvnames{g}) = [];
        end
        max_qtildel = [];
        fweights = [];
    end
    
    % Retrieve maximum lag order for serial correlation test
    psc_du2 = xlsread(intfname,'DOMINANT UNIT','psc_du_2');    
    
    
    if gxvnum == 1 % univariate model
        % retrieve the estimation type: in levels or in first differences
        [junk esttype_du_char] = xlsread(intfname,'DOMINANT UNIT','esttype_du_2'); %#ok
        if strcmp(esttype_du_char,'levels');
            esttype_du2 = 0;
        elseif strcmp(esttype_du_char,'first differences');
            esttype_du2 = 1;
        end              
        
        % retrieve the case of treatment of deterministics 
        estcase_du2 = xlsread(intfname,'DOMINANT UNIT','deterministics_du_2');

        % note that if the model is estimated in first differences, only
        % the intercept can be included, ensure this is the case below
        if esttype_du2 == 1
            estcase_du2 = 0;
        end
    else
        esttype_du2 = [];
        estcase_du2 = [];
    end
    
   
    
    if isfeedback == 0 && gxvnum > 1
        % no augmented regression, recover the VAR coefficients of the
        % (multivariate) dominant unit model from the VECM estimates
        
        % Set the lag order of the augmented GVAR
        %*****************************************
        amaxlag = max([maxlag varxlag_du]);
        
        % recover the VAR coefficients 
        [a0_du a1_du Theta_du] =  vec2var_du(amaxlag,gxvnum,varxlag_du,alpha_du,beta_du,Psi_du,estcase_du);
        maxfeednum = length(vnames);
        Lambda_du = zeros(gxvnum,maxfeednum,amaxlag);
        
        ptildel_found = [];
        qtildel_found = [];
    else
        % Estimation of the augmented regression
        %***************************************
        if lagselect_du2 == 1
            
            if pauseflag == 1
                % Lag order selection
                [ptildel_found qtildel_found tab] = augmentedregression(isfeedback,lagselect_du2,info2,...
                    max_ptildel,max_qtildel,maxlag,gxv,gxvnames,exog_du,exognames_du,feedbacks_flagmat,ecm_du,alpha_du,beta_du,estcase_du,esttype_du2,estcase_du2,psc_du2);
                
                
                % print to interface file the resulting lag orders
                xlswrite(intfname,ptildel_found,'DOMINANT UNIT','Q5');
                if isfeedback == 1
                    for i=1:length(isfeedback_byeq)
                        if isfeedback_byeq(i) == 0
                            qtildel_found(i) = NaN; 
                        end
                    end
                    xlswrite(intfname,qtildel_found,'DOMINANT UNIT','AE5');
                end
                
                if printout == 1
                    % print results of lag order selection in output.xlsx
                    disp('- Adding to output.xlsx: Lag order selection results for the dominant unit model augmented regression');
                    disp('  equation(s) and residual serial correlation test results');
                    xlswrite([outdir 'output.xlsx'],tab,'DUmodel_augres_sc');
                end
                
                % pause and let user eventually modify the lag orders for the
                % augmented regression
                %************************************************
                winopen(intfname);
                disp(' ');
                msg = sprintf('>>> Pause and go to %s: Check the lag orders found for the augmented regression',intfname);
                disp(msg);
                disp('     (or inputted if no selection criterion was previously used), then press enter.');
                pause
                disp(' ');
                msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                disp(msg);
                endpause = input(' ','s');
                disp('The program is now running (do not press any key)  ');
                disp(' ');
            end
            
            % retrieve lag orders from cells found_ptildel found_qtildel
            ptildel_found = xlsread(intfname,'DOMINANT UNIT','found_ptildel');      
            qtildel_found_tmp = xlsread(intfname,'DOMINANT UNIT','found_qtildel');
            qtildel_found = qtildel_found_tmp(:,1);
            if length(qtildel_found) < gxvnum
                qtildel_found = [qtildel_found; NaN(gxvnum-length(qtildel_found),1)];
            end
            
        else
            ptildel_found = max_ptildel;
            if isfeedback == 1 
                qtildel_found = max_qtildel;
            else
                qtildel_found = [];
            end
        end
        
        % Estimation of the augmented regression
        [junk junk tab a0_du a1_du Theta_du Lambda_du amaxlag] = augmentedregression(isfeedback,0,[],...
            ptildel_found,qtildel_found,maxlag,gxv,gxvnames,exog_du,exognames_du,feedbacks_flagmat,ecm_du,alpha_du,beta_du,estcase_du,esttype_du2,estcase_du2,psc_du2); %#ok
        
        %
        
        if printout == 1
            % print OLS estimates of augmented regression in output.xlsx
            disp('- Adding to output.xlsx: Estimates of the dominant unit model (Augmented regression)');
            xlswrite([outdir 'output.xlsx'],tab,'DUmodel_augres_est');
            
            if isfeedback == 1
            % print the dominant unit model data, including xtilde
            % variables
                disp('- Adding to countrydata.xlsx: Data for the dominant unit model including the feedback variables');
            else
                disp('- Adding to countrydata.xlsx: Data for the dominant unit model');
            end
            vlab = date;
            hlab = {'date'};
            tab = [];
            for i=1:length(gxvnames)
                hlab = [hlab gxvnames(i)]; %#ok
                tab = [tab gxv(:,i)]; %#ok
            end
            
            if isfeedback == 1
                vartilde = zeros(1,length(vnames));
                xtildeblock = [];
                xtildelabel = {};
                vflag = sum(feedbacks_flagmat,1);
                for j=1:length(vnames)
                    if vflag(j) > 0
                        varname = vnames{j};
                        for i=1:rows(feedbacks_flagmat)
                            if feedbacks_flagmat(i,j) == 1
                                ind = i;
                                varpos = 0;
                                for k = 1:length(exognames_du.(gxvnames{ind}))
                                    if strcmp(exognames_du.(gxvnames{ind})(k),varname)
                                        varpos = k;
                                    end
                                end
                                break
                            end
                        end
                        xtildeblock = [xtildeblock exog_du.(gxvnames{ind})(:,varpos)]; %#ok
                        out = sprintf('%stilde',varname);
                        out = {out};
                        xtildelabel = [xtildelabel out]; %#ok
                    end
                end
            end
            
            if isfeedback == 0
                toprint = [hlab; vlab num2cell(tab)];
            elseif isfeedback == 1
                toprint = [[hlab xtildelabel]; vlab num2cell(tab) num2cell(xtildeblock)];
            end
            
            xlswrite([outdir 'countrydata.xlsx'],toprint,'DOMINANT UNIT');
        end
    end
    
  
    
   
    % Create country-specific lists of global exogenous variables which enter as weak
    % exogenous variables (useful for computing the augmented link matrices and persistence profiles)
    for i = 1:cnum
        gxpointer.(cnames{i}) = zeros(1,gxvnum);
        gxcount.(cnames{i}) = [];
        for j = 1:gvnum
            if gvflag(i,j) == 1 % global variable j is included in country i
                % is global variable j a global exogenous variable?
                for k=1:length(gxvnames)
                    if strcmp(gvnames{j},gxvnames{k})
                        % yes, add it in counter
                        gxpointer.(cnames{i})(k) = 1;
                    end
                end
            end
        end
        gxcount.(cnames{i}) = sum(gxpointer.(cnames{i}));
    end    
    
    % include list of global exogenous variables in the endoglist structure
    endoglist.du_model = gxvnames;
    
    % print empty cell on MAIN sheet so to ensure that next time Matlab
    % opens the interface file, it opens at the MAIN sheet
    xlswrite(intfname,NaN,'MAIN','C1');
    
else
    a0_du = [];
    a1_du = [];
    Theta_du = [];
    Lambda_du = [];
    maxlag_du = [];
    varxlag_du = [];
    rank_du = [];
    beta_du = [];
    beta_norm_du = [];
    beta_norm.du_model = [];
    estcase_du = [];
    esttype_du2 = [];
    estcase_du2 = [];
    psc_du2 = [];
    ptildel_found = [];
    qtildel_found = [];
    gxidx = [];
    isfeedback = [];
    feedbacks_flagmat = [];
    fweights = [];
    exog_du = [];
    exognames_du = [];
    gxpointer = [];
    gxcount = [];
    amaxlag = maxlag;
end




%% 4) SOLVING THE GVAR MODEL
disp(' ');
disp('4) SOLVING THE GVAR MODEL');
disp('********************************************************************');

%% 4.1) Creating link matrices W_i
%**************************************************************************
disp('4.1) Creating the link matrices W_i')

% Creating W_i link matrices, z_it vectors and the  x_t global vector
%**********************************************************************


% retrieve the specific year of the weights matrix at which the gvar has to
% be solved. Clearly, if weights are fixed, the year is irrelevant
if strcmp(weightstype_flag,'fixed')
    yearpos = 1;
elseif strcmp(weightstype_flag,'time-varying')
    tvw_solution = xlsread(intfname,'MAIN','tvw_solution');
    yearpos = find(yearseq_estimation==tvw_solution);
end

% define the matrix which is gonna be used to solve the GVAR
wmat_sol = wm_t.(ylabelseq_estimation{yearpos});

% create link matrices
[z x xnames W] = create_linkmatrices(cnames,cnum,vnum,vnames,fvnum, fvnames, gnvnum,gnvnames,endog,endoglist,exognx,exognxlist,dvtype,wmat_sol);



% create the y vector, which includes the x vector and the global
% exogenous variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = [x; gxv'];
ynames = [xnames; gxvnames'];

Ky = rows(y); % number of vars in the GVAR (including possible dominant units)
K = Ky - gxvnum; % number of endogenous variables in the GVAR
 


if not(isempty(gxv))
    % retrieve also the z vector augmented by the global (exogenous)
    % variables, denote it zy
    zy = create_linkmatrices(cnames,cnum,vnum,vnames,fvnum, fvnames, gnvnum,gnvnames,endog,endoglist,exog,exoglist,dvtype,wmat_sol);
    
    % augment the W matrices accordingly
    for n=1:cnum
        Wy.(cnames{n}) = [W.(cnames{n}) zeros(rows(W.(cnames{n})),gxvnum); zeros(gxvnum,K) diag(gxpointer.(cnames{n}))];
    end
    
    % create also an augmented W matrix for the dominant unit model
    Wy.du_model = [zeros(gxvnum,K) eye(gxvnum)];
else
    zy = z;
    Wy = W;
end




% Creating country numbers
%*************************
cnumbers = zeros(K,1);

pos=1;
for n=1:cnum

    tmp = length(endoglist.(cnames{n}));
    cnumbers(pos:pos+tmp-1) = n;

    pos = pos+tmp;
end

% Creating a K x 1 cell of country names
%***************************************
cnames_x = {};
cnames_s_x = {};

for i=1:K
    cnames_x{i} = cnames_long{cnumbers(i)}; %#ok
    cnames_s_x{i} = cnames{cnumbers(i)}; %#ok
end

cnames_x = cnames_x';
cnames_s_x = cnames_s_x';

cnames_y = cnames_x;
cnames_s_y = cnames_s_x;
for i=K+1:Ky
    cnames_y{i} = 'DOMINANT UNIT MODEL';
    cnames_s_y{i} = 'du_model';
end


if dumodel_flag == 1 && isfeedback == 1
    % Create the matrix Wtilda of feedback variables of the dominant unit model
    %**************************************************************************
    maxfeednum = max(sum(feedbacks_flagmat,2));
    Wtilde = zeros(maxfeednum,K);
    k=0;
    for g=1:cols(feedbacks_flagmat)
        if sum(feedbacks_flagmat(:,g)) > 0
            k=k+1;
            for j=1:K
                if isfield(fweights.(vnames{g}),cnames_s_x{j}) && strcmp(xnames{j},vnames{g})
                    Wtilde(k,j) = fweights.(vnames{g}).(cnames_s_x{j});
                end
            end
        end
    end
else
    Wtilde = [];
end



if forc_flag == 1 && misal == 1     
    % creating also a global vector with remaining actual data for comparing it with GVAR
    % forecasts
    rx=[];
    for n=1:cnum
        endblock_tmp = rendog.(cnames{n});
        endblock = endblock_tmp';
        rx = [rx; endblock]; %#ok
    end
    
    if dumodel_flag == 1 % add also the remaining data for the dominant unit global variables
        for i=1:gxvnum
            rx = [rx; rgv.(gxvnames{i})']; %#ok
        end
    end
    
    clear endblock endblock_tmp
end




%% 4.2) Retrieving the parameters of the VARX* models from the VECMX* estimates
%**************************************************************************
disp('4.2) Retrieving the parameters of the VARX* models from the VECMX* estimates');


[a0 a1 Theta Lambda0_tmp Lambda_tmp] =  vecx2varx(maxlag,cnum,cnames,zy,endoglist,varxlag,alpha,beta,Psi,estcase);

lsize = size(Lambda_tmp.(cnames{1}));
if length(lsize) == 3
    nlag = lsize(3);
elseif length(lsize) == 2
    nlag = 1;
end

for n=1:cnum
    Lambda0.(cnames{n}) = [];
    Lambda.(cnames{n}) = [];
  
    Gamma0.(cnames{n}) = zeros(length(endoglist.(cnames{n})),gxvnum);
    Gamma.(cnames{n}) = zeros(length(endoglist.(cnames{n})),gxvnum,nlag);    
    
    for j=1:length(exoglist.(cnames{n}))
        flag = 0;
        for i=1:length(gxvnames)
            if strcmp(exoglist.(cnames{n})(j),gxvnames{i})
                flag = 1;    
                Gamma0.(cnames{n})(:,i) = Lambda0_tmp.(cnames{n})(:,j);
                Gamma.(cnames{n})(:,i,:) = Lambda_tmp.(cnames{n})(:,j,:);
            end
        end
        if flag == 0;
            Lambda0.(cnames{n}) = [Lambda0.(cnames{n}) Lambda0_tmp.(cnames{n})(:,j)];
            Lambda.(cnames{n}) = [Lambda.(cnames{n}) Lambda_tmp.(cnames{n})(:,j,:)];
        end
    end
end

%% 4.3) Constructing and solving the GVAR model
%**************************************************************************
disp('4.3) Constructing and solving the GVAR model ');


[delta_0 delta_1 H0 C zeta Sigma_zeta eta Sigma_eta eigens mods mlag] = solve_GVAR(maxlag,amaxlag,...
                          cnum,cnames,W,a0,a1,Theta,Lambda0,Lambda,Gamma0,Gamma,x,gxv,Wtilde,a0_du,a1_du,...
                          Theta_du,Lambda_du); %#ok
% note: in the code above, we are correctly using the W matrices (and not
% Wy, which are augmented for global exog variables)                          
out = sprintf('- Lag order of the GVAR: %d', mlag); disp(out);


if printout == 1
    disp('- Adding to output.xlsx: Eigenvalues of the GVAR');
    % print eigenvalues and moduli
    
    title = {'Eigenvalues of the GVAR Model in Descending Order'};
    cotitle = {'Corresponding Moduli'};
    tab = [title num2cell(NaN(1,5)) cotitle];
    tab = [tab; num2cell(NaN(1,7))];
    
    eig_out = {};
    realp = [];
    imagp = [];

    for i=1:length(eigens)
        realp = [realp; real(eigens(i))]; %#ok
        imagp = [imagp; imag(eigens(i))]; %#ok
    end

    [realp_sorted idx] = sort(realp,'descend');

    imagp_sorted = imagp(idx);

    for i=1:length(eigens)
        if imagp_sorted(i) == 0
            eig_idv = sprintf('%-3.14f', realp_sorted(i));
        else
            eig_idv = sprintf('%-3.14f %+3.14fi', realp_sorted(i), imagp_sorted(i));
        end
        eig_out = [eig_out; eig_idv]; %#ok
    end

    tab = [tab; eig_out num2cell(NaN(rows(eig_out),5)) num2cell(mods)];
    
    xlswrite([outdir 'output.xlsx'],tab,'eigenval');
    
    clear realp realp_sorted
end


% cleaning
clear title cotitle eig_out mods imagp imagp_sorted fv gv eigens cnumbers

disp(' ');
disp('End of the GVAR model estimation.');
disp(' ');

%% 4.4) GVAR forecasts
%**************************************************************************

%================================================================
Rmin=0.25; % lower bound to be imposed on interest rate (ex-ante) forecasts  
trIR=(1/freq)*log(1+(Rmin/100)); %change here if a different tranformation is
                                 %required for the interest rate variables 
%================================================================


if forc_flag == 1

    disp('4.4) GVAR forecasts '); 
 
    % retrieving interest rate variables for imposing lower bounds 
     
     [junk irlb] = xlsread(intfname,'MAIN','vforlb');  
     
     if isempty(irlb)==0;
     irlbstr=irlb{1};
     virlb=convertstr2vec(irlbstr);
     else
     virlb=irlb;    
     end
     
     % Additional lower bound restrictions on forecasts
    lb_flag = xlsread(intfname,'MAIN','lb_flag');
    
    fhorz = xlsread(intfname,'MAIN','fhorz');

 if isempty(irlb)==0

    %% Lower bound restrictions on the ex-ante forecasts of the interest rates
    %**************************************************************************
    disp('- Imposing lower bound restrictions on the ex-ante forecasts of the interest rates');

 end

if lb_flag==1

    %% Impose (additional) lower bound restrictions on the ex-ante forecasts 
    %**************************************************************************
    disp('- Imposing (additional) lower bound restrictions on the ex-ante forecasts');
    
    lb_restr = ExF_lb_restr(cnames_y, ynames,outdir);
else
   lb_restr=[]; 
end
    
    
     yforc = forecast_GVAR(fhorz,sposfirst,sposlast,Ky,y,mlag,C,delta_0,delta_1,...
         trIR,virlb,ynames,lb_flag,lb_restr);


    if printout == 1
        disp('- Adding to output.xlsx: GVAR (ex-ante) forecasts');
  
        if annual == 0
            misal = not(strcmp(lastobs,max_date));
        else
            misal = not(isequal(lastobs,max_date));
        end
        
        if not(misal)
            rx = [];
        end
        
  
        if cols(rx) > fhorz % if horizon of actual data is larger than forecast horizon, trim actual data
            print_forecasts(lastobs,fhorz,cnames_s_y,ynames,yforc,rx(:,1:fhorz),outdir);
        else
            print_forecasts(lastobs,fhorz,cnames_s_y,ynames,yforc,rx,outdir);
        end
    end
end


if con_forc_flag==1
    disp('4.4b) GVAR conditional forecasts ');
    
    dummy_check = 1;
    while dummy_check==1
         con_fhorz = xlsread(intfname,'MAIN','con_fhorz'); % H
         con_fhorz_restr= xlsread(intfname,'MAIN','con_fhorz_restr'); % H_bar
        
        if con_fhorz_restr < con_fhorz    
            disp('');
            disp('Warning: The restriction horizon must be greater or equal to the forecast horizon.');
            disp('');

            % pause and let user change the restriction horizon
            %***************************************************
            winopen(intfname);
            disp(' ');
            msg = sprintf('>>>> Pause and go to %s: Change the restriction horizon so that it is greater  ',intfname);
            disp(msg);
            disp('        or equal to the forecast horizon.');
            pause
            disp(' ');
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s');
            disp('The program is now running (do not press any key)  ');
            disp(' ');

        elseif  con_fhorz_restr >= con_fhorz
            dummy_check = 0;
        end 
    end 
    
   con_fdate = getForecastDateSeries(lastobs,con_fhorz_restr); 
   
  disp('- Imposing the conditional forecast restrictions');
   [con_forc_restr_mtx con_forc_restr_x] = con_forc_restr(Ky,con_fhorz_restr,ynames,lastobs,cnames_s_y,outdir) ;
   
   [mu_s] = con_forecast_GVAR(mlag,sposfirst,sposlast,...
    delta_0,delta_1,y,Ky,Sigma_eta,C,...
    con_fhorz,con_fhorz_restr,...
    con_forc_restr_mtx,con_forc_restr_x,ynames);


    if printout == 1
        disp('- Adding to output.xlsx: GVAR conditional forecasts');

        
        if annual == 0
            misal = not(strcmp(lastobs,max_date));
        else
            misal = not(isequal(lastobs,max_date));
        end
        
        if not(misal)
            rx = [];
        else
            % tebaldi: adding else because rx is not defined
            
            rx = [];
        end
        
        if cols(rx) > con_fhorz % if horizon of actual data is larger than forecast horizon, trim actual data
            print_con_forecasts(lastobs,con_fhorz,cnames_s_y,ynames,mu_s,rx(:,1:con_fhorz),outdir);
        else
            print_con_forecasts(lastobs,con_fhorz,cnames_s_y,ynames,mu_s,rx,outdir);
        end
    end 
end 
    
  
  
%% 4.5) Trend/Cycle (TC) decomposition of the GVAR model
%**************************************************************************
if TCflag == 1
    disp('4.5) Trend/Cycle(TC) decomposition of the GVAR model');

% redefining and making sure that the esttype and estcase are not empty (required for the TC decomposition) 
  if isempty(estcase_du)==1
  estcase_du_TC=-9999;
  else
  estcase_du_TC=estcase_du;    
  end   

  if isempty(esttype_du2)==1
  esttype_du2_TC=-9999; 
  else
  esttype_du2_TC=esttype_du2;    
  end  

  if isempty(estcase_du2)==1
  estcase_du2_TC=-9999;
  else
  estcase_du2_TC=estcase_du2;    
  end
    
  TC_RestrictionFlag = xlsread(intfname, 'MAIN','TC_RestrictionFlag');% Flag for whether to put trend restrictions

if dumodel_flag==0
%Case 2 for all countries  
cond_TC_noDU=sum(estcase_tmp)==2*size(estcase_tmp,1);   
else
cond_TC_noDU=-9999;
end

if dumodel_flag==1
    %Case 2 for VEC augmented model or univariate model in levels with no trend
    cond_TC_DU=sum(estcase_tmp)==2*size(estcase_tmp,1) && (estcase_du_TC==2 || (esttype_du2_TC==0 && estcase_du2_TC==0));
else
    cond_TC_DU=-9999;
end

if TC_RestrictionFlag==1
   if (dumodel_flag==0 && cond_TC_noDU==1)||(dumodel_flag==1 && cond_TC_DU==1);
        % pause and make user change the trend restriction flag
        %***************************************************
        winopen(intfname);
        disp(' ');
        msg = sprintf('>>>> Pause and go to %s: Based on your current settings no trends are',intfname);
        disp(msg);
        disp('        included in the GVAR model and the trend restrictions in the T/C decomposition ');
        disp('        field should be disabled (i.e. set to zero). ');
        pause
        disp(' ');
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
        
        TC_RestrictionFlag = xlsread(intfname, 'MAIN','TC_RestrictionFlag');% Flag for whether to put trend restrictions
    end
end

if TC_RestrictionFlag==1
    [TC_restr]= TC_trend_restr(cnames_y,ynames,outdir) ;
    while sum(TC_restr)==0;
       disp(' ');
       disp(' >>> Warning: Based on the selected settings, trend restrictions need to be imposed for the T/C decomposition.');
       disp('     Further instructions are given below.');
       disp(' ');   
    [TC_restr]= TC_trend_restr(cnames_y,ynames,outdir) ;
        if sum(TC_restr)~=0;
        break;
        end   
    end        
else
    TC_restr = false(length(cnames_y),1);
end

[yp yc y_tilde]= TCdecomp(TC_RestrictionFlag,C,zeta,mlag,y,TC_restr,cond_TC_noDU,cond_TC_DU,dumodel_flag);


disp('- Writing to TCdecomposition.xls: Domestic variables, permanent components and cyclical');
disp('  components (deviations from steady states)');
print_TCdecomp(y, yp, yc, y_tilde, vnames,gvnames,ynames, date, ...
    cnames_y,mlag, outdir);
end
clear TC_restr yp yc y_tilde yp_dt yp_st a
%**************************************************************************


if strcmp(intfname,'gvarBriefDemo.xls')
    do_da = 1;
else
    disp(' ');
    disp(' ');
    disp('Do you want to go on with the dynamic analysis of the GVAR? If yes, type y, otherwise type n:');
    disp(' ');
    disp('======================================================================================');
    disp('Dynamic analysis of the GVAR includes: Computation of the persistence profiles, impulse');
    disp('response analysis, forecast error variance decomposition, and computation of the critical');
    disp('values for the structural stability tests and the logLik test for overidentifying '); 
    disp('restrictions on the cointegrating vectors. The last two assuming that the corresponding');
    disp('functions, as well as the bootstrap, are enabled.');
    disp('======================================================================================');
    disp('');
    
    
    k = 1;
    while k < 2
        answer_tmp = input(' ','s');
        answer = {answer_tmp};
        
        if strcmp(answer,'y')
            do_da = 1;
            k=2;
        elseif strcmp(answer,'n')
            do_da = 0;
            k=2;
        end
    end
end

%% **************************************************************************
%**************************************************************************


if do_da == 1


    pack storevars % (technical, manages stored objects for memory purposes)


    %% 5) IMPULSE RESPONSE ANALYSIS

    disp(' ');
    disp('5) IMPULSE RESPONSE ANALYSIS');
    disp('********************************************************************');


    %% 5.1) Settings
    %************************************************************************************************
 disp('5.1) Settings, shock selection and (optional) region definition for regional GIRFs and GFEVDs');
    
    
    % check whether a dominant unit model exists
    if dumodel_flag == 1
        cnamesy_long = [cnames_long; 'DOMINANT UNIT MODEL'];
        cnamesy = [cnames; 'du_model'];
    else
        cnamesy_long = cnames_long;
        cnamesy = cnames;
    end     
 
    % check whether any specification is inputted
    vshockmtrx = xlsread(intfname,'MAIN','vshockmtrx');
    
    if isempty(vshockmtrx) % no shocks specified, force pauseflag =1
        pauseflag = 1;
    end
    

    if pauseflag == 1

        % retrieve past info about data availability for each model
        % for GIRFs & GFEVDs part

        dvflag_da = zeros(rows(dvflag),cols(dvflag));
        if not(gvnum==0)
            gvflag_da = zeros(rows(gvflag),cols(gvflag));
        end

        for i=1:cnum
            for j=1:cols(dvflag)
                if dvflag(i,j) == 1
                    dvflag_da(i,j) = 0;
                elseif dvflag(i,j) == 0;
                    dvflag_da(i,j) = NaN;
                elseif isnan(dvflag(i,j))
                    dvflag_da(i,j) = NaN;
                end
            end

            if not(gvnum==0)
                for g=1:cols(gvflag)
                    if gvflag(i,g) == 1
                        gvflag_da(i,g) = NaN;
                    elseif gvflag(i,g) == 0
                        gvflag_da(i,g) = NaN;
                    elseif isnan(gvflag(i,g))
                        gvflag_da(i,g) = NaN;
                    elseif gvflag(i,g) == 2
                        gvflag_da(i,g) = 0;
                    end
                end
            end
        end
   
        if not(gvnum==0) && not(isempty(gxv))
            % add cells for the global exogenous variables
            gvflag_da = [gvflag_da; NaN(1,gvnum)];
            for j=1:gvnum
                if gxidx(j) == 1
                    gvflag_da(end,j) = 0;
                end
            end
        end
        
        if not(gvnum==0)
            xlswrite(intfname,gvflag_da,'MAIN','FO5');
        end
        xlswrite(intfname,dvflag_da,'MAIN','ET5');

        xlswrite(intfname,[cnamesy_long cnamesy],'MAIN','EQ5');
        
        % pause and let user specify settings
        %**************************************************************************
        winopen(intfname);
        disp(' ');
        msg = sprintf('>>> Pause and go to %s. Define:',intfname);
        disp(msg);   
        disp(' ');
        disp('- The settings for the dynamic analysis');
        disp('- The country and associated variable ordering, if structural GIRFs (SGIRFs)');
        disp('  or orthogonalised IRFs (OIRFs) and corresponding FEVDs are to be carried out');
        disp('- The shocks to be performed');
        disp('- The region composition, if you wish to carry out regional shocks and/or obtain');
        disp('  regional impulse response functions and forecast decomposition results. ');
        disp(' ');
               disp('  Once the above are defined, press enter.');
        disp(' ');
        disp('======================================================================================= ');
        disp('Input at this stage is required, unless you are using the brief demo interface file.');
        disp('1. For selecting the covariance matrix for computation of the point (and bootstrap)')
        disp('   estimates make sure you enter 1 in one of the three fields of your choice, and');
        disp('   that the other two fields are set to 0.');
        disp('2. Under the bootstrap approach "inverse", and the "shuffle" approach when combined with');
        disp('   the selection of OIRFs for the impulse response functions, the program will allow you');
        disp('   to change your selections for performing or not shrinkage on the correlation matrix');
        disp('   as it continues to run. This will also be the case under OIRFs when no bootstrap is');
        disp('   enabled.');
        disp('3. Under the bootstrap approach "shuffle" when combined with GIRFs or SGIRFs, information');
        disp('   in the shrinkage panel is not applicable and any selected options will be ignored. By'); 
        disp('   default no shrinkage will be performed.');
        disp('4. If you are computing structural GIRFs which permit structural identification of shocks'); 
        disp('   to a SINGLE country, in the adjacent panel define the country of interest using the');
        disp('   country short name and specify the ordering of the endogenous variables of that country'); 
        disp('   (only those that are included in the GVAR model). Do the same for all countries if you');
        disp('   are computing orthogonalised IRFs. In both cases if a dominant unit is present, this');
        disp('   should be placed first and its variables also ordered (see the user guide for details).') 
        disp('5. Define the shock(s) you wish to entertain. You can select country, regional and/or');
        disp('   global shocks. Click on the "Select shocks" heading in the interface file for further');
        disp('   guidance. Please refer to the user guide for more details.');
        disp('6. If you wish to entertain regional shocks and/or to obtain regional impulse responses');
        disp('   and forecast error decompositions using the PPP-GDP data provided by the program (or');
        disp('   any other data provided for aggregation), scroll to the far right of the MAIN worksheet');
        disp('   and define your desired regions (if not already defined) by consulting the user guide.');
        disp('7. For every run of the program the shocks will need to be redefined. To retain the');
        disp('   defined shocks in a subsequent run, disable the "Run the program with pauses" function');
        disp('   at the initial settings stage.');
        disp(' ');
        disp(' USING THE FULL DEMO INTERFACE FILE');
        disp('If you are using the full demo interface file and would like to replicate the results');
        disp('in the Output Full Demo folder, define:                   ');
        disp(' 1. A regional positive shock to real GDP in the "Rest of Asia" region (enter 2 in the');
        disp('    real GDP cell of, for example, Korea or any other country that belongs to that region).');
        disp(' 2. A global negative shock to real equity prices (enter -3 in the equity cell of, for');
        disp('    example, Australia, or in any other country).');
        disp(' 3. A positive shock to US short-term interest rate (enter 1 in the short-term interest');
        disp('    rate cell of the US model).');
        disp(' 4. A positive shock to the oil price (enter 1 in the oil price cell of the dominant unit');
        disp('    model.');
        disp('======================================================================================== ');
        disp(' ');
        pause
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
    end

    % Settings
    %*********
    
    bootstrap_flag = xlsread(intfname,'MAIN','bootstrap_flag'); % =1 if bootstrap is performed

    if bootstrap_flag == 1
        B = xlsread(intfname,'MAIN','B'); % # of bootstrap samples to generate
        
        %==================================================================
        vtrunc = floor(0.05*B);
        while vtrunc==0
            disp(' ');
            disp(' >>> Warning: The number of bootstrap replications specified in the dynamic analysis section');
            disp('     of the MAIN sheet of the interface file needs to be at least 20. Please make the appropriate');
            disp('     adjustment.');
           
            disp(' ');

            winopen(intfname);
            pause
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s'); 
            disp('The program is now running (do not press any key)  ');
            disp(' ');
            B = xlsread(intfname,'MAIN','B');
            vtrunc=floor(0.05*B);
        end
        %==================================================================
    end

    N = xlsread(intfname,'MAIN','N'); % horizon's length for PPs, GIRFs and GFEVDs
    
    [junk sgirfflag_name] = xlsread(intfname,'MAIN','sgirfflag'); %#ok
    % Set whether you want to do Generalized IRFs, Structural Generalized IRFs
    % or Orthogonalized IRFs
    if strcmp('GIRFs',sgirfflag_name) 
        sgirfflag=0;  
    elseif strcmp('SGIRFs',sgirfflag_name) 
        sgirfflag=1;   
    elseif strcmp('OIRFs',sgirfflag_name) 
        sgirfflag=2;         
    end
    

    
    
    %% 5.2) Calculating varcov matrices
    %**************************************************************************
    disp('5.2) Computing the covariance matrix');

    % retrieve user choices:

    % for point estimates (the same covariance structure is used for
    % bootstrap estimates and for generating the bootstrap data)
    pe_samplecov =  xlsread(intfname,'MAIN','pe_samplecov');
    pe_blockdiag =  xlsread(intfname,'MAIN','pe_blockdiag');
    pe_blockdiagexc =  xlsread(intfname,'MAIN','pe_blockdiagexc');
    
    vcheck = (pe_samplecov==0) && (pe_blockdiag==0) && (pe_blockdiagexc==0);
    

    if vcheck ~= 0
        % pause and let user specify settings
        %**************************************************************************
        winopen(intfname);
        disp(' ');
        msg = sprintf('>>> Warning: Make sure you have specified the type of covariance matrix');
        disp(msg);
        disp('    you wish to use for the computation of the point estimates and the bootstrap');
        disp('    (if enabled), by entering 1 in one of the three fields of your choice and 0')
        disp('    in the remaining two fields, and then press enter');
        disp(' ');
        pause
        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
        disp(msg);
        endpause = input(' ','s');
        disp('The program is now running (do not press any key)  ');
        disp(' ');
    end
        

    % Transformation of the varcov matrix 
    %**********************************************************************    
    pe_samplecov =  xlsread(intfname,'MAIN','pe_samplecov');
    pe_blockdiag =  xlsread(intfname,'MAIN','pe_blockdiag');
    pe_blockdiagexc =  xlsread(intfname,'MAIN','pe_blockdiagexc');
    [junk pe_country_exc] =  xlsread(intfname,'MAIN','pe_exccountry'); %#ok
    
    if pe_samplecov == 1
        % use original estimated varcov matrix
        pe_meth = 1;
    elseif pe_blockdiag == 1
        % transform the original varcov matrix to a block diagonal one
        pe_meth = 2;
    elseif pe_blockdiagexc == 1
        % transform the original varcov matrix to a block diagonal one,
        % while leaving unrestricted the cross-country correlations for a
        % given exception country
        pe_meth = 3;
    end    
    
    
    % Settings for shrinkage
    %**********************************************************************
    use_shrinkedvcv =  xlsread(intfname,'MAIN','use_shrinkedvcv');
    use_shrinkedvcv_dg =  xlsread(intfname,'MAIN','use_shrinkedvcv_dg');

    if bootstrap_flag== 1
        % Recover info about the type of bootstrap
        [junk Bootstrap_approach] = xlsread(intfname,'MAIN','Bootstrap_approach'); %#ok
        
        if strcmp('inverse',Bootstrap_approach) == 1
            shuffleflag=0;
        elseif strcmp('shuffle',Bootstrap_approach) == 1
            shuffleflag=1;
        end
    else
        shuffleflag = [];
    end
       
    

    
    % Transform the varcov matrix 
    pe_varcov_tx = transform_varcov(pe_meth,pe_country_exc,Sigma_zeta,cnames_s_y); 
    
    
    
    % Shrinkage of variance covariance matrix and check whether resulting varcov is positive
    % definite (if doing boostrap or orthogonalized IRFs)
    %**************************************************************
    
    arbitrarylambda = []; do_alambda = 0;
    if pauseflag == 0 % no pauses
        % check if the user has inputed an arbitrary value for the
        % shrinkage parameter lambda
        arbitrarylambda = xlsread(intfname,'MAIN','arbitrarylambda');
        
        if not(isempty(arbitrarylambda))
            do_alambda = 1;
        end
    end
    
    if bootstrap_flag == 0
        if sgirfflag == 0 || sgirfflag == 1 % GIRFs, SGIRFs
            if use_shrinkedvcv == 0 % don't do shrinkage
                pe_varcov=pe_varcov_tx;
            elseif use_shrinkedvcv == 1 % do shrinkage
                if do_alambda == 1
                    lambda_star = arbitrarylambda;
                    pe_varcov =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                else
                    [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                    
                    % print the resulting shrinkage parameter to the
                    % interface file
                    xlswrite(intfname,lambda_star,'MAIN','EH33');
                    
                end
            end
        elseif sgirfflag == 2 % OIRFs
            % check if the pe_varcov_tx matrix is positive definite
            posflag=1;
            try
                chol(pe_varcov_tx);
            catch %#ok
                posflag=0;
            end
            firstpause = 1;
            
            if posflag == 1
                if use_shrinkedvcv == 0
                    % don't do shrinkage
                    pe_varcov=pe_varcov_tx;
                elseif use_shrinkedvcv == 1
                    % shrinkage
                    if do_alambda == 1
                        lambda_star = arbitrarylambda;
                        pe_varcov =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                    else
                        [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                        
                        % print the resulting shrinkage parameter to the
                        % interface file
                        xlswrite(intfname,lambda_star,'MAIN','EH33');
                        
                    end
                end
            elseif posflag == 0
                k=1;
                while k<2
                    
                    
                    disp('>>> Warning: The covariance matrix is not positive definite.')
                    disp(' ');
                    if firstpause == 1 || (firstpause == 0 && isempty(arbitrarylambda))
                        msg = sprintf('>>> Pause and go to %s: Ensure that perform shrinkage on the correlation matrix',intfname);
                        disp(msg);
                        disp('         for computing the point and bootstrap estimates is set to 1 (if it is not already),');
                        disp('         then press enter.');
                    elseif firstpause == 0 && not(isempty(arbitrarylambda))
                        msg = sprintf('>>> Pause and go to %s: Increase the value of the shrinkage parameter',intfname);
                        disp(msg);
                        disp('         without making any other changes, then press enter.');
                    end
                    disp(' ');
                    winopen(intfname);
                    pause
                    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                    disp(msg);
                    disp('The program is now running (do not press any key)  ');
                    disp(' ');
                    
                    arbitrarylambda =  xlsread(intfname,'MAIN','arbitrarylambda');
                    use_shrinkedvcv =  xlsread(intfname,'MAIN','use_shrinkedvcv');
                    firstpause = 0;
                    
                    
                    if use_shrinkedvcv == 1
                        if do_alambda == 1
                            lambda_star = arbitrarylambda;
                            pe_varcov =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                        else
                            [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                            
                            % print the resulting shrinkage parameter to the
                            % interface file
                            xlswrite(intfname,lambda_star,'MAIN','EH33');
                            
                            k = 2; % as with internally computed lambda, always varcov is positive definite
                        end
                    else
                        pe_varcov=pe_varcov_tx;
                    end
                    
                    posflag=1;
                    try
                        chol(pe_varcov);
                    catch %#ok
                        posflag=0;
                    end
                    
                    if posflag == 1
                        k = 2;
                        disp('- OK: The covariance matrix is positive definite');
                    end
                end
            end
        end
    elseif bootstrap_flag == 1
        if shuffleflag == 1 % shuffle-option of bootstrap
            
            if sgirfflag == 0 || sgirfflag == 1 % GIRFs, SGIRFs
                % no shrinkage
                pe_varcov = pe_varcov_tx;
                pe_varcov_dg = pe_varcov_tx;
                lambda_star = [];
                
            elseif sgirfflag == 2 % OIRFs
                % check if the pe_varcov_tx matrix is positive definite
                posflag=1;
                try
                    chol(pe_varcov_tx);
                catch %#ok
                    posflag=0;
                end
                firstpause = 1;
                
                if posflag == 1
                    % don't do shrinkage
                    pe_varcov = pe_varcov_tx;
                    pe_varcov_dg = pe_varcov_tx;
                    lambda_star = [];
                    
                elseif posflag == 0
                    k=1;
                    while k<2
                        
                        disp('>>> Warning: The covariance matrix is not positive definite.')
                        disp(' ');
                        if firstpause == 1 || (firstpause == 0 && isempty(arbitrarylambda))
                            msg = sprintf('>>> Pause and go to %s: Ensure that perform shrinkage on the correlation matrix',intfname);
                            disp(msg);
                            disp('        for generating the bootstrap data and for computing the point and bootstrap');
                            disp('        estimates are set to 1 (if they are not already), then press enter.');
                        elseif firstpause == 0 && not(isempty(arbitrarylambda))
                            msg = sprintf('>>> Pause and go to %s: Increase the value of the shrinkage parameter',intfname);
                            disp(msg);
                            disp('         without making any other changes, then press enter.');
                        end
                        disp(' ');
                        winopen(intfname);
                        pause
                        msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                        disp(msg);
                        disp('The program is now running (do not press any key)  ');
                        disp(' ');
                        
                        arbitrarylambda =  xlsread(intfname,'MAIN','arbitrarylambda');
                        use_shrinkedvcv =  xlsread(intfname,'MAIN','use_shrinkedvcv');
                        use_shrinkedvcv_dg =  xlsread(intfname,'MAIN','use_shrinkedvcv_dg');
                        firstpause = 0;
                        
                        
                        if use_shrinkedvcv == 1
                            if do_alambda == 1
                                lambda_star = arbitrarylambda;
                                pe_varcov =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                            else
                                [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                                
                                % print the resulting shrinkage parameter to the
                                % interface file
                                xlswrite(intfname,lambda_star,'MAIN','EH33');
                                
                                k = 2; % as with internally computed lambda, always varcov is positive definite
                            end
                        elseif use_shrinkedvcv == 0
                            % the varcov matrix for point and bootstrap estimates is not positive definite and
                            % the user does not want to use shrinkage
                            pe_varcov=pe_varcov_tx;
                        end
                        
                        if use_shrinkedvcv_dg == 1
                            if do_alambda == 1
                                pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                            else
                                [pe_varcov_dg lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                                
                                % print the resulting shrinkage parameter to the
                                % interface file
                                xlswrite(intfname,lambda_star,'MAIN','EH33');
                                
                                k = 2;
                            end
                        elseif use_shrinkedvcv_dg == 0
                            pe_varcov_dg = pe_varcov_tx;
                        end
                        
                        posflag=1;
                        if use_shrinkedvcv == 1
                            try
                                chol(pe_varcov);
                            catch %#ok
                                posflag=0;
                            end
                        end
                        
                        if posflag == 1
                            k = 2;
                            if use_shrinkedvcv == 1
                                disp('- OK: The covariance matrix is positive definite');
                            end
                        end
                    end
                end
            end
        elseif shuffleflag == 0 % inverse-option of bootstrap, shrinkage
            % check if the pe_varcov_tx matrix is positive definite
            posflag=1;
            try
                chol(pe_varcov_tx);
            catch %#ok
                posflag=0;
            end
            firstpause = 1;
            
            if posflag == 1
                if use_shrinkedvcv == 0
                    % don't do shrinkage
                    pe_varcov = pe_varcov_tx;
                elseif use_shrinkedvcv == 1
                    % shrinkage
                    if do_alambda == 1
                        lambda_star = arbitrarylambda;
                        pe_varcov =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                    else
                        [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                        
                        % print the resulting shrinkage parameter to the
                        % interface file
                        xlswrite(intfname,lambda_star,'MAIN','EH33');

                    end
                end
                % note: the shrinkage of the varcov matrix for bootstrap data
                % generation is dealt after the loop
                
            elseif posflag == 0
                k=1;
                while k<2
                    
                    disp('>>> Warning: The covariance matrix is not positive definite.')
                    disp(' ');
                    if firstpause == 1 || (firstpause == 0 && isempty(arbitrarylambda))
                        msg = sprintf('>>> Pause and go to %s: Ensure that perform shrinkage on the correlation matrix',intfname);
                        disp(msg);
                        disp('         for generating the bootstrap data is set to 1 (if it is not already),');
                        disp('         then press enter.');
                        disp(' ');
                        disp('       - If GIRFs or SGIRFs are selected, at this stage you can choose to retain or');
                        disp('         change your selection associated with performing shrinkage on the correlation');
                        disp('         matrix for computing point and bootstrap estimates.');
                        disp(' ');
                        disp('       - If OIRFs are selected, ensure in addition that perform shrinkage on the');
                        disp('         correlation matrix for point and bootstrap estimates is set to 1 (if it is not already).');
                    elseif (firstpause == 0 && not(isempty(arbitrarylambda)))
                        msg = sprintf('>>> Pause and go to %s: Increase the value of the shrinkage parameter',intfname);
                        disp(msg);
                        disp('         without making any other changes, then press enter.');
                    end
                    disp(' ');
                    winopen(intfname);
                    pause
                    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                    disp(msg);
                    disp('The program is now running (do not press any key)  ');
                    disp(' ');
                    
                    arbitrarylambda =  xlsread(intfname,'MAIN','arbitrarylambda');
                    use_shrinkedvcv =  xlsread(intfname,'MAIN','use_shrinkedvcv');
                    use_shrinkedvcv_dg =  xlsread(intfname,'MAIN','use_shrinkedvcv_dg');
                    firstpause = 0;
                    
                    if use_shrinkedvcv == 0 && use_shrinkedvcv_dg == 1 % possible for GIRFs or SGIRFs
                        % the varcov matrix for point and bootstrap estimates is not positive definite and
                        % the user does not want to use shrinkage
                        pe_varcov = pe_varcov_tx;
                        
                        % but the varcov matrix for bootstrap data
                        % generation gets shrinked
                        if do_alambda == 1
                            lambda_star = arbitrarylambda;
                            pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                        else
                            [pe_varcov_dg lambda_star] =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                            
                            % print the resulting shrinkage parameter to the
                            % interface file
                            xlswrite(intfname,lambda_star,'MAIN','EH33');
                            
                            k = 2;
                        end
                    elseif use_shrinkedvcv == 1 && use_shrinkedvcv_dg == 1
                        if do_alambda == 1
                            lambda_star = arbitrarylambda;
                            pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                            pe_varcov = pe_varcov_dg;
                        else
                            [pe_varcov_dg lambda_star] =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                            pe_varcov = pe_varcov_dg;
                            
                            % print the resulting shrinkage parameter to the
                            % interface file
                            xlswrite(intfname,lambda_star,'MAIN','EH33');
                            
                            k = 2;
                        end
                    end
                    
                    posflag=1;
                    try
                        chol(pe_varcov_dg);
                    catch %#ok
                        posflag=0;
                    end
                    
                    if posflag == 1
                        k = 2;
                        if use_shrinkedvcv == 1 || use_shrinkedvcv_dg == 1
                            disp('- OK: The covariance matrix is positive definite');
                        end
                    end
                end
            end
        end
    end
    
    
    % select resulting varcov matrices
    if use_shrinkedvcv == 0
        pe_varcov = pe_varcov_tx;
        lambda_star = [];
    end
    
    if bootstrap_flag == 1 && shuffleflag == 0 % note that the case of shuffleflag == 1 is dealt above
        if use_shrinkedvcv_dg == 0 % don't shrink varcov matrix for bootstrap data generation
            pe_varcov_dg = pe_varcov_tx;
        elseif use_shrinkedvcv_dg == 1
            if do_alambda == 0
                [pe_varcov_dg lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                
                % print the resulting shrinkage parameter to the
                % interface file
                xlswrite(intfname,lambda_star,'MAIN','EH33');
                
            elseif do_alambda == 1 && pauseflag == 0
                lambda_star = arbitrarylambda;
                pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
            end
        end
    end
    
    
    if printout == 1;
        % print information about shrinkage of varcov matrix
        tab = ['Shrinkage of Covariance Matrix' num2cell(NaN(1,1))];
        tab = [tab; num2cell(nan(1,2))];
        tab = [tab; 'Covariance matrix for point and bootstrap estimates:' num2cell(NaN(1,1))];
        if use_shrinkedvcv == 0 || isempty(lambda_star);
            tab = [tab; '- used original selected covariance matrix (no shrinkage performed)' num2cell(NaN(1,1))];
        elseif use_shrinkedvcv == 1
            tab = [tab; '-- shrinkage was performed' num2cell(NaN(1,1))];
            if do_alambda == 0
                tab = [tab; '-- internally computed shrinkage parameter =' num2cell(lambda_star)];
            elseif do_alambda == 1
                tab = [tab; '-- shrinkage parameter provided by the user =' num2cell(lambda_star)];
            end
        end
        if bootstrap_flag == 1 && shuffleflag == 0
            tab = [tab; num2cell(nan(1,2))];
            tab = [tab; 'Covariance matrix for bootstrap data generation:' num2cell(NaN(1,1))];
            if use_shrinkedvcv_dg == 0 || isempty(lambda_star);
                tab = [tab; '- used original selected covariance matrix (no shrinkage performed)' num2cell(NaN(1,1))];
            elseif use_shrinkedvcv_dg == 1
                tab = [tab; '-- shrinkage was performed' num2cell(NaN(1,1))];
                if do_alambda == 0
                    tab = [tab; '-- internally computed shrinkage parameter =' num2cell(lambda_star)];
                elseif do_alambda == 1
                    tab = [tab; '-- shrinkage parameter provided by the user =' num2cell(lambda_star)];
                end
            end
        end
        xlswrite([outdir 'output.xlsx'],tab,'cov_shrinkage');
    end
    

    
    
    
    
    if sgirfflag == 1 || sgirfflag == 2 % you want to do Structural GIRF or OIRFs
        % note: this part of code handles both cases of SGIRFs and OIRFs
        
        % retrieve info about first countries
        %************************************
        
        [junk firstcountries] = xlsread(intfname,'MAIN','firstcountries'); %#ok
        [junk newordervars] = xlsread(intfname,'MAIN','newordervars'); %#ok
        
        firstcountries = firstcountries(not(strcmp(firstcountries,'')));
        
        [pe_varcov_s Sigma_zeta0 H0_s_t C_s Wy_s beta_s beta_norm_s ncnames ncnames_long ...
            nendoglist nynames ncnames_y ncnames_s_y yorder nyorder] = reorder_GVAR(firstcountries,...
         newordervars,cnamesy,cnamesy_long,Ky,endoglist,...
         H0,C,mlag,ynames,cnames_s_y,Wy,beta,beta_norm,estcase,pe_varcov);
        
       if dumodel_flag == 1
           ncnamesy = ['du_model'; ncnames];
           ncnamesy_long = ['DOMINANT UNIT MODEL'; ncnames_long];
       else
           ncnamesy = ncnames;
           ncnamesy_long = ncnames_long;
       end    
    else
        Sigma_zeta0 = [];
        H0_s_t = [];
        firstcountries = [];
        newordervars = [];
        yorder = [];
        nyorder = [];
    end
    
    
    if sgirfflag == 1 || sgirfflag == 2
        % save original matrices of countries and variables original ordering
        cnames_backup = cnames;
        cnames_long_backup = cnames_long;
        cnamesy_backup = cnamesy;
        cnamesy_long_backup = cnamesy_long;
        endoglist_backup = endoglist;
        ynames_backup = ynames;
        cnames_y_backup = cnames_y;
        cnames_s_y_backup = cnames_s_y;
    end    
    
    
    %% 5.3) Calculation of country and regional PPP-GDP weights
    %**************************************************************************
    disp('5.3) Calculation of country and regional PPP-GDP weights');

    cweights = country_weights(aggrwgts,cnames,vnum,vnames,gnvnum,gnvnames,endoglist,1);

    if printout == 1;
        % print country weights
        %**********************
        disp('- Adding to output.xlsx: Country weights');
        print_cweights(cnames,cnames_long,vnames,gnvnames,cweights,outdir);
    end


    % check whether there exist regions for aggregation of IRFs and FEVDs
    %****************************************************************
    [junk rnames_long] = xlsread(intfname,'MAIN','rnames_long'); %#ok % region names, format long
    [junk rnames] = xlsread(intfname,'MAIN','rnames'); % region names, format short

    if isempty(rnames_long)
        % user does not specify regions for aggregation of GIRFs
        regres_flag = 0;
        regions = [];
        rweights = [];
    else
        % specify regions
        regres_flag = 1;

        rnames_long = rnames_long(not(strcmp(rnames_long,'')));
        rnames = rnames(not(strcmp(rnames,'')));

        [junk regcountries] = xlsread(intfname,'MAIN','regcountries'); %#ok


        regcountries_temp = [];
        k=1;
        for i=1:length(regcountries)
            if strcmp(regcountries{i},'')
                regions.(rnames{k}) = regcountries_temp;
                k=k+1;
                regcountries_temp = [];
            else
                regcountries_temp = [regcountries_temp regcountries(i)]; %#ok
                if k==length(rnames)
                    regions.(rnames{k}) = regcountries_temp;
                end
            end
        end

        regcountries = regcountries(not(strcmp(regcountries,'')));

        clear junk regcountries_temp
    end
    
    if regres_flag == 1
        % Calculation of regional weights
        %*********************************
        disp('- Adding to output.xlsx: Regional weights');
        rweights = regional_weights(vnames, gnvnames, rnames, regions, cweights);

        if printout == 1
            % print regional weights
            %***********************
            print_rweights(vnames,gnvnames,rnames,regions,cweights,rweights,outdir);
        end
    end


    %% 5.4) Point estimates
    %**************************************************************************
    disp('5.4) Point estimates');

    
    % Calculating dynamic multipliers
    %********************************
    disp('- Computing the dynamic multipliers');
    if sgirfflag == 1 || sgirfflag == 2
        % compute two sets of matrices, one reordered (for SGIRFs/OIRFs),
        % the other not reordered (for Persistence Profiles)
        PHI = dyn_multipliers(Ky,mlag,C,N); % not reordered
        PHI_s = dyn_multipliers(Ky,mlag,C_s,N); % reordered
    else
        PHI = dyn_multipliers(Ky,mlag,C,N); % not reordered
    end
        

    % 5.4a) Persistence Profiles
    %**************************************************************************
    disp('- Computing the persistence profiles');
    
    
    if isempty(beta_norm.du_model)
        cnamesy_pp = cnamesy(not(strcmp(cnamesy,'du_model')));
        cnamesy_long_pp = cnamesy_long(not(strcmp(cnamesy_long,'DOMINANT UNIT MODEL')));
    else
        cnamesy_pp = cnamesy;
        cnamesy_long_pp = cnamesy_long;
    end
    
    if dumodel_flag == 1
        % augment beta_norm matrices with zeros for excluded global exogenous
        % variables, so that matrices are conformable in computing persistence
        % profiles
        for n=1:length(cnamesy)
            if not(strcmp(cnamesy(n),'du_model'))
                % by construction the dominant unit model contains all global
                % exogenous variables, so skip it
                oldblock = beta_norm.(cnamesy{n})(end+1-gxcount.(cnamesy{n}):end,:);
                newblock = zeros(gxvnum,rank.(cnamesy{n}));
                kk = 1;
                for k=1:gxvnum
                    if gxpointer.(cnamesy{n})(k) == 1
                        newblock(k,:) = oldblock(kk,:);
                        kk = kk + 1;
                    end
                end
                beta_norm_gx.(cnamesy{n}) = [beta_norm.(cnamesy{n})(1:end-gxcount.(cnamesy{n}),:); newblock];
            end
        end
        beta_norm_gx.du_model = beta_norm.du_model;
    else
        beta_norm_gx = beta_norm;
    end
    
    
    
    PPres_t = pprofile(PHI,pe_varcov,H0,Wy,beta_norm_gx,N,cnamesy_pp,estcase); % for persistence profiles, use not reordered matrices
    
    PPres = PPres_t';

    k=1;
    for n=1:length(cnamesy_pp)
        PP.(cnamesy_pp{n}) = [];
        for j=1:rank.(cnamesy_pp{n})
            PP.(cnamesy_pp{n}) = [PP.(cnamesy_pp{n}) PPres(:,k)];
            k=k+1;
        end
    end

    
    if printout == 1
        
        disp('- Adding to output.xlsx: Persistence profiles');
        % print PProfiles
        title = {'Persistence Profile of the Effect of System-Wide Shocks to the Cointegrating Relations of the GVAR Model'};

        hlabel1 = [];
        hlabel2 = [];
        table = [];
        for n=1:length(cnamesy_pp)
            for j=1:cols(PP.(cnamesy_pp{n}))
                hlabel1 = [hlabel1 cnamesy_long_pp(n)]; %#ok
                tmp = sprintf('CV%d',j);
                tmp = {tmp};
                hlabel2 = [hlabel2 tmp]; %#ok
                table = [table PP.(cnamesy_pp{n})(:,j)]; %#ok
            end
        end

        vlabel1 = {'Horizon'};
        vlabel2 = (0:N)';
        
        tab = [title num2cell(NaN(1,cols(table)))];
        tab = [tab; num2cell(NaN(1,1+cols(table)))];
        tab = [tab; num2cell(NaN(1,1)) hlabel1];
        tab = [tab; vlabel1 hlabel2];
        tab = [tab; num2cell([vlabel2 table])];
        
        xlswrite([outdir 'output.xlsx'],tab,'PP');
           
    end
    


    % 5.4b) IRFs and associated FEVDs
    %**************************************************************************
    disp('- IRFs and FEVDs');
   
   
    vshockmtrx = xlsread(intfname,'MAIN','vshockmtrx');
    if not(isempty(gxv)) % add row of NaNs for the dominant unit model
        vshockmtrx = [vshockmtrx; NaN(1,vnum)];
    end
    
    if not(gvnum==0)

        % use a trick, as in matlab 2010: xlsread ignores any outer rows or
        % columns of the spreadsheet that contain no numeric data. If there
        % are single or multiple nonnumeric rows at the top or bottom, or
        % single or multiple nonnumeric columns to the left or right,
        % xlsread does not include these rows or columns in the output.

        gvshockmtrx_tmp_a = xlsread(intfname,'MAIN','gvshockmtrx');
        gvshockmtrx_tmp_b = gvshockmtrx_tmp_a(:,2:end);
        gvshockmtrx_tmp = [gvshockmtrx_tmp_b; NaN(cnum-rows(gvshockmtrx_tmp_b),cols(gvshockmtrx_tmp_b))];
        
        gvshockmtrx = zeros(rows(vshockmtrx),gvnum);
        
        for j=1:gvnum
            for i=1:rows(gvshockmtrx)
                if j <= cols(gvshockmtrx_tmp)
                    if not(isnan(gvshockmtrx_tmp(i,j))) && not(gvshockmtrx_tmp(i,j) ==0)
                        % 18/02/2012: fixed problem of sign of global shock
                        % gvshockmtrx(i,j) = 1;
                        gvshockmtrx(i,j) = gvshockmtrx_tmp(i,j);
                    end
                end
            end
        end
    end
    
    
    if not(gvnum==0)
        shockmtrx = [vshockmtrx(:,1:vnum) gvshockmtrx(:,1:gvnum)];
        allvnames = [vnames gvnames];
    else
        shockmtrx = vshockmtrx(:,1:vnum);
        allvnames = vnames;
    end


    %%%%%
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
                    countype = [];
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
                    countype = [];
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
                


                
                if sgirfflag == 1 || sgirfflag == 2 % use reordered matrices PHI_s, pe_varcov_s, H0_s_t
                    
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
                    
                    IRFMAT = irf(Ky,N,PHI_s,pe_varcov_s,H0_s_t,eslct,sgirfflag,Sigma_zeta0);
                    FEVDMAT = fevd(Ky,N,PHI_s,pe_varcov_s,H0_s_t,eslct,sgirfflag,Sigma_zeta0);           
                else
                    IRFMAT = irf(Ky,N,PHI,pe_varcov,H0,eslct,sgirfflag,Sigma_zeta0);
                    FEVDMAT = fevd(Ky,N,PHI,pe_varcov,H0,eslct,sgirfflag,Sigma_zeta0);
                end
                
                if printout == 1

                    simtitle = sprintf('- Shock simulation no.%d: Writing results to output folder',ind_sim);
                    disp(simtitle);
                    
                    if sgirfflag == 1 || sgirfflag == 2
                        % use new country and variable ordering for
                        % printing output
                        cnames = ncnames;
                        cnames_long = ncnames_long;
                        cnamesy = ncnamesy;
                        cnamesy_long = ncnamesy_long;
                        endoglist = nendoglist;
                        ynames = nynames;
                        cnames_y = ncnames_y;
                        cnames_s_y = ncnames_s_y;
                    end
                    
                    
                    disp('- Creating the output file irfs.xls');
                    % printing IRFs
                    irfdir = print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,IRFMAT, ...
                        shocktype,shocksign,vartype,countype,regtype,0,0,sgirfflag,'Point Estimates');
                    movefile('Output\irfs.xls',irfdir);

                    if graphsflag == 1
                        disp('- Plotting IRFs');
                        % plotting IRFs
                        plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                            irfdir,shocktype,shocksign,N,Ky,cnames_y,ynames,IRFMAT,'graphs','Point Estimates',1,0,sgirfflag);
                    end

                    disp('- Creating the output file fevds.xls');
                    % printing FEVDs
                    fevddir = print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,FEVDMAT, ...
                        vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,...
                        vartype,countype,regtype,outdir_t,0,0,sgirfflag,'Point Estimates');
                    movefile('Output\fevds.xls',fevddir);

                    if regres_flag == 1 % if Region-specific

                        % aggregate country responses into regional
                        % responses
                        rIRFMAT = aggregate_results(rweights,vnames,gnvnames,rnames,regions,ynames,cnames_s_y,IRFMAT);
                        % print regional responses
                        disp('- Creating the output file rirfs.xls');
                        print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,rIRFMAT, ...
                            shocktype,shocksign,vartype,countype,regtype,0,1,sgirfflag,'Point Estimates');
                        movefile('Output\rirfs.xls',irfdir);

                        if graphsflag == 1
                            plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                irfdir,shocktype,shocksign,N,Ky,cnames_y,ynames,rIRFMAT,'graphs','Point Estimates',1,1,sgirfflag);
                        end

                        % aggregate country FEVDs into regional
                        rFEVDMAT = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,FEVDMAT);
                        % printing regional FEVDs
                        disp('- Creating the output file rfevds.xls');
                        fevddir = print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,rFEVDMAT, ...
                            vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,...
                            vartype,countype,regtype,outdir_t,0,1,sgirfflag,'Point Estimates');
                        movefile('Output\rfevds.xls',fevddir);
                    end
                    
                    if sgirfflag == 1 || sgirfflag == 2
                        % set back original ordering of countries and
                        % variables
                        cnames = cnames_backup;
                        cnames_long = cnames_long_backup;
                        cnamesy = cnamesy_backup;
                        cnamesy_long = cnamesy_long_backup;
                        endoglist = endoglist_backup;
                        ynames = ynames_backup;
                        cnames_y = cnames_y_backup;
                        cnames_s_y = cnames_s_y_backup;
                    end 
                end
            end
        end
    end



    if bootstrap_flag == 1

        % 5.5) GVAR Bootstrap
        %**************************************************************************
        disp('5.5) GVAR Bootstrap');
        

        % Bootstrap:  
      [median_PP lbound_PP ubound_PP median_IRF lbound_IRF ubound_IRF median_FEVD lbound_FEVD ubound_FEVD ind_sim overid_LR_95cv overid_LR_99cv] = bootstrap_GVAR(B,N, ...
           zeta,pe_meth,lambda_star,use_shrinkedvcv,pe_varcov_dg,pe_country_exc,cnames_s_y,y,ynames,cnamesy,cnamesy_long,cnum,cnames,cweights,rnames,regions,rweights,...
           delta_0,delta_1,C,H0,K,dvtype,ntypes,wmatrices,maxlag,amaxlag,mlag,varxlag,estcase,rank,beta_r,...
          dv, vnames, vnames_long, fvnames, gvnames, gvnames_long, vnum, fvnum, gvnum, gnvnum, gnvnames, gxvnames, gxidx,gxcount,gxpointer, maxlag_du, varxlag_du, rank_du, estcase_du,...
          isfeedback,feedbacks_flagmat,fweights,exog_du,exognames_du,esttype_du2,estcase_du2,psc_du2,ptildel_found,qtildel_found,Wtilde,...
          dvflag, fvflag, gvflag, allvnames, shockmtrx,...
          ydate,ylabelseq_estimation,yearseq_estimation,nyears_estimation,nobs,wmat_sol,redoflag,sgirfflag,firstcountries,newordervars,overid_flag,shuffleflag,yorder,nyorder);

       
      % transpose pprofiles
      median_PPt = median_PP';
      lbound_PPt = lbound_PP';
      ubound_PPt = ubound_PP';
      
      k=1;
      for n=1:length(cnamesy_pp)
          PPmed.(cnamesy_pp{n}) = [];
          PPlb.(cnamesy_pp{n}) = [];
          PPub.(cnamesy_pp{n}) = [];
          for j=1:rank.(cnamesy_pp{n})
              PPmed.(cnamesy_pp{n}) = [PPmed.(cnamesy_pp{n}) median_PPt(:,k)];
              PPlb.(cnamesy_pp{n}) = [PPlb.(cnamesy_pp{n}) lbound_PPt(:,k)];
              PPub.(cnamesy_pp{n}) = [PPub.(cnamesy_pp{n}) ubound_PPt(:,k)];
              k=k+1;
          end
      end
      
        if printout == 1
           
            disp('- Adding to output.xlsx: Median estimates and confidence bands for persistence profiles');
            % print PProfiles

            title1 = {'Persistence Profile of the Effect of System-Wide Shocks to the Cointegrating Relations of the GVAR Model - Bootstrap Median estimates'};
            title2 = {'Persistence Profile of the Effect of System-Wide Shocks to the Cointegrating Relations of the GVAR Model - Bootstrap Lower bounds'};
            title3 = {'Persistence Profile of the Effect of System-Wide Shocks to the Cointegrating Relations of the GVAR Model - Bootstrap Upper bounds'};

            hlabel1 = [];
            hlabel2 = [];
            tab1 = [];
            tab2 = [];
            tab3 = [];
            for n=1:length(cnamesy_pp)
                for j=1:cols(PP.(cnamesy_pp{n}))
                    hlabel1 = [hlabel1 cnamesy_long_pp(n)]; %#ok
                    tmp = sprintf('CV%d',j);
                    tmp = {tmp};
                    hlabel2 = [hlabel2 tmp]; %#ok
                    tab1 = [tab1 PPmed.(cnamesy_pp{n})(:,j)]; %#ok
                    tab2 = [tab2 PPlb.(cnamesy_pp{n})(:,j)]; %#ok
                    tab3 = [tab3 PPub.(cnamesy_pp{n})(:,j)]; %#ok
                end
            end

            vlabel1 = {'Horizon'};
            vlabel2 = (0:N)';

            % print median estimates
            tab = [title1 num2cell(NaN(1,cols(tab1))); num2cell(NaN(1,1+cols(tab1)))];
            tab = [tab; num2cell(NaN(1,1)) hlabel1; vlabel1 hlabel2; num2cell([vlabel2 tab1])];
            xlswrite([outdir 'output.xlsx'],tab,'PP_bs_median');

            % print lower bounds
            tab = [title2 num2cell(NaN(1,cols(tab2))); num2cell(NaN(1,1+cols(tab2)))];
            tab = [tab; num2cell(NaN(1,1)) hlabel1; vlabel1 hlabel2; num2cell([vlabel2 tab2])];            
            xlswrite([outdir 'output.xlsx'],tab,'PP_bs_lbound');

            % print upper bounds
            tab = [title3 num2cell(NaN(1,cols(tab3))); num2cell(NaN(1,1+cols(tab3)))];
            tab = [tab; num2cell(NaN(1,1)) hlabel1; vlabel1 hlabel2; num2cell([vlabel2 tab3])];              
            xlswrite([outdir 'output.xlsx'],tab,'PP_bs_ubound');
            
        end


        if printout == 1

            if overid_flag == 1
                disp('- Adding to output.xlsx: Overidentifying restrictions test statistic and bootstrapped critical values');
                % printing overidentifying restrictions tests statistic and
                % critical values
                %*************************************************************
                title = {'Test of Overidentifying Restrictions: Statistics and Bootstrap Critical Values'};
                
                hlabel = {'Country' 'Likelihood Ratio statistic' 'Degrees of freedom' '95% critical values' '99% critical values'} ;
                vlabel = cnames_long;
                
                table = [];
                for n=1:cnum
                    table = [table; overid_LR.(cnames{n}) overid_dgf.(cnames{n}) overid_LR_95cv.(cnames{n}) overid_LR_99cv.(cnames{n})]; %#ok
                end

                tab = [title num2cell(NaN(1,cols(table)))];
                tab = [tab; num2cell(NaN(1,1+cols(table)))];
                tab = [tab; hlabel];
                tab = [tab; vlabel num2cell(table)];
                
                xlswrite([outdir 'output.xlsx'],tab,'overid_restr_logLik');
                
                clear title hlabel vlabel table tab
            end

            

            % printing IRFs and FEVDs
            %**************************************************************************
            disp('- Writing IRF and FEVD results to output folder');

            ind_sim=0;
            for j=1:cols(shockmtrx)
                for n=1:rows(shockmtrx)
                    if (not(shockmtrx(n,j) == 0)) && (not(isnan(shockmtrx(n,j))))
                        ind_sim = ind_sim+1;
                        simtitle = sprintf('- Shock simulation no.%d: Writing results to output folder',ind_sim);
                        disp(simtitle);

                        shocksign = sign(shockmtrx(n,j));
                        shocktype = abs(shockmtrx(n,j));


                        if shocktype == 1  % country-specific shock
                            regtype = [];
                            countype = cnamesy(n);
                            vartype = allvnames(j);

                        elseif shocktype == 2
                            countype = [];
                            country = cnamesy(n);
                            for r=1:length(rnames)
                                for l=1:length(regions.(rnames{r}))
                                    if strcmp(country,regions.(rnames{r})(l))
                                        regtype = rnames(r);
                                    end
                                end
                            end
                            vartype = allvnames(j);


                        elseif shocktype == 3  % global shock
                            regtype=[];
                            countype = [];
                            vartype = allvnames(j);
                            countrylist = fieldnames(cweights.(allvnames{j}));
                        end


                        % printing IRFs and FEVDs
                        %*************************
                        if sgirfflag == 1 || sgirfflag == 2
                            % use new country and variable ordering for
                            % printing output
                            cnames = ncnames;
                            cnames_long = ncnames_long;
                            cnamesy = ncnamesy;
                            cnamesy_long = ncnamesy_long;
                            endoglist = nendoglist;
                            ynames = nynames;
                            cnames_y = ncnames_y;
                            cnames_s_y = ncnames_s_y;
                        end
                        
                        disp('- Creating the output file irfs_bs.xls');
                        irfdir_bs = print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,median_IRF(:,:,ind_sim), ...
                            shocktype,shocksign,vartype,countype,regtype,1,0,sgirfflag,'Median Estimates');
                        print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,lbound_IRF(:,:,ind_sim), ...
                            shocktype,shocksign,vartype,countype,regtype,1,0,sgirfflag,'Lower Bounds');
                        print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,ubound_IRF(:,:,ind_sim), ...
                            shocktype,shocksign,vartype,countype,regtype,1,0,sgirfflag,'Upper Bounds');
                        movefile('Output\irfs_bs.xls',irfdir_bs);

                        if graphsflag == 1
                            disp('- Plotting IRFs');
                            plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,median_IRF(:,:,ind_sim),'graphs_bs','Median Estimates',1,0,sgirfflag);
                            plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,lbound_IRF(:,:,ind_sim),'graphs_bs','Lower Bounds',0,0,sgirfflag);
                            plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,ubound_IRF(:,:,ind_sim),'graphs_bs','Upper Bounds',0,0,sgirfflag);
                        end

                        disp('- Creating the output file fevds_bs.xls');
                        fevddir_bs = print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,median_FEVD(:,:,ind_sim), ...
                            vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,0,sgirfflag,'Median Estimates');
                        print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,lbound_FEVD(:,:,ind_sim), ...
                            vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,0,sgirfflag,'Lower Bounds');
                        print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,ubound_FEVD(:,:,ind_sim), ...
                            vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,0,sgirfflag,'Upper Bounds');
                        movefile('Output\fevds_bs.xls',fevddir_bs);

                        % IRFs can be country-specific and region-specific

                        if regres_flag == 1 % if Region-specific

                            % aggregate country responses into regional
                            % responses
                            median_rIRF = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,median_IRF(:,:,ind_sim));
                            lbound_rIRF = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,lbound_IRF(:,:,ind_sim));
                            ubound_rIRF = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,ubound_IRF(:,:,ind_sim));

                            % print regional responses
                            disp('- Creating the output file rirfs_bs.xls');
                            irfdir_bs = print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,median_rIRF, ...
                                shocktype,shocksign,vartype,countype,regtype,1,1,sgirfflag,'Median Estimates');
                            print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,lbound_rIRF, ...
                                shocktype,shocksign,vartype,countype,regtype,1,1,sgirfflag,'Lower Bounds');
                            print_irfs(outdir_t,cnamesy,cnamesy_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,ynames,cnames_y,N,ubound_rIRF, ...
                                shocktype,shocksign,vartype,countype,regtype,1,1,sgirfflag,'Upper Bounds');
                            movefile('Output\rirfs_bs.xls',irfdir_bs);

                            if graphsflag == 1
                                plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                    irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,median_rIRF,'graphs_bs','Median Estimates',1,1,sgirfflag);
                                plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                    irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,lbound_rIRF,'graphs_bs','Lower Bounds',0,1,sgirfflag);
                                plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnamesy,countype,cnamesy_long,...
                                    irfdir_bs,shocktype,shocksign,N,Ky,cnames_y,ynames,ubound_rIRF,'graphs_bs','Upper Bounds',0,1,sgirfflag);
                            end

                            % aggregate country-level FEVDs into regional
                            % responses %
                            median_rFEVD = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,median_FEVD(:,:,ind_sim));
                            lbound_rFEVD = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,lbound_FEVD(:,:,ind_sim));
                            ubound_rFEVD = aggregate_results(rweights,vnames,gvnames,rnames,regions,ynames,cnames_s_y,ubound_FEVD(:,:,ind_sim));

                            % print regional FEVDs
                            disp('- Creating the output file rfevds_bs.xls');
                            fevddir_bs = print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,median_rFEVD, ...
                                vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,1,sgirfflag,'Median Estimates');
                            print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,lbound_rFEVD, ...
                                vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,1,sgirfflag,'Lower Bounds');
                            print_fevds(cnamesy,cnamesy_long,rnames,rnames_long,ynames,cnames_y,N,ubound_rFEVD, ...
                                vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,vartype,countype,regtype,outdir_t,1,1,sgirfflag,'Upper Bounds');
                            movefile('Output\rfevds_bs.xls',fevddir_bs);
                        end
                        
                        if sgirfflag == 1 || sgirfflag == 2
                            % set back original ordering of countries and
                            % variables
                            cnames = cnames_backup;
                            cnames_long = cnames_long_backup;
                            cnamesy = cnamesy_backup;
                            cnamesy_long = cnamesy_long_backup;
                            endoglist = endoglist_backup;
                            ynames = ynames_backup;
                            cnames_y = cnames_y_backup;
                            cnames_s_y = cnames_s_y_backup;
                        end
                    end
                end
            end   
        end
    end

    disp('--- End of the dynamic analysis ---');
end


if ss_flag == 1
    %% Bootstrap critical values for structural stability tests
    ss_bootstrap_flag = xlsread(intfname,'MAIN','ss_bootstrap_flag');
    
    if ss_bootstrap_flag == 1
        disp(' ');
        disp('5.6) Bootstrap critical values for structural stability tests');
        
        ss_B = xlsread(intfname,'MAIN','ss_B');
        
        vtrunc = int16(0.05*ss_B);
        while vtrunc==0
            disp(' ');
            disp(' >>> Warning: The number of bootstrap replications for computation of the critical values');
            disp('     of the structural stability tests needs to be at least 20. Please make the appropriate');
            disp('     adjustment within the settings section of the MAIN sheet of the interface file.');
            disp(' ');
            
            winopen(intfname);
            pause
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s');
            disp('The program is now running (do not press any key)  ');
            disp(' ');
            ss_B = xlsread(intfname,'MAIN','ss_B');
            vtrunc=int16(0.05*ss_B);
        end
        %==================================================================
        
        
        
        
        if pauseflag == 1
            % pause and let user specify settings
            %**************************************************************************
            winopen(intfname);
            disp(' ');
            msg = sprintf('>>> Pause and go to %s: Select the required bootstrap approach, the covariance',intfname);
            disp(msg);
            disp('     matrix to be used for generating the bootstrap data and whether to perform shrinkage on');
            disp('     the corresponding correlation matrix, and then press enter.');
            disp(' ');
            disp('    (If these are already defined and you do not wish to change them, simply close');
            disp('     the interface file and press enter).');
            disp(' ');
            disp('======================================================================================= ');
            disp(' Notes:                                                                                   ');
            disp('i.   For the selection of the covariance matrix choose one of the given three matrices ');
            disp('     by selecting 1 in the field of your choice, and making sure that the other fields ');
            disp('     are set to 0. ');
            disp('ii.  The field related to shrinkage performance for the point and bootstrap estimates is');
            disp('     not applicable when bootstraping to obtain the critical values for the structural');
            disp('     stability tests. Any selected options will be ignored.');
            disp('iii. Under the bootstap approach "inverse" the program will allow you to change your ');
            disp('     selection for performing or not shrinkage on the correlation matrix used to generate');
            disp('     the bootstrap data, when it continues to run.');
            disp('IV.  Under the bootstrap approach "shuffle", information in the shrinkage panel is not');
            disp('     applicable and any selected options will be ignored. By default no shrinkage');
            disp('     will be performed.');
            disp(' ');
            disp(' USING THE FULL DEMO INTERFACE FILE');
            disp('If you are using the full demo interface file and would like to replicate the results');
            disp('in the Output Full Demo folder simply retain the selected settings.                  ');
            disp('======================================================================================== ');
            disp(' ');
            pause
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s');
            disp('The program is now running (do not press any key)  ');
            disp(' ');
        end
        
        disp('- Computing the covariance matrix');
        
        % retrieve user choice for bootstrap varcov matrix:
        
        pe_samplecov =  xlsread(intfname,'MAIN','pe_samplecov');
        pe_blockdiag =  xlsread(intfname,'MAIN','pe_blockdiag');
        pe_blockdiagexc =  xlsread(intfname,'MAIN','pe_blockdiagexc');
        [junk pe_country_exc] =  xlsread(intfname,'MAIN','pe_exccountry'); %#ok
        
        
        vcheck = (pe_samplecov==0) && (pe_blockdiag==0) && (pe_blockdiagexc==0);
        if vcheck ~= 0
            % pause and let user specify settings
            %**************************************************************************
            winopen(intfname);
            disp(' ');
            msg = sprintf('>>> Warning: Make sure you have specified the type of covariance matrix');
            disp(msg);
            disp('    you wish to use for performing the bootstrap by entering 1 in one of the three')
            disp('    fields of your choice, and 0 in the remaining two fields, and then press enter.');
            
            disp(' ');
            pause
            msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
            disp(msg);
            endpause = input(' ','s');
            disp('The program is now running (do not press any key)  ');
            disp(' ');
        end
        
        
        % Transformation of the varcov matrix for bootstrap
        %******************************************************************
        pe_samplecov =  xlsread(intfname,'MAIN','pe_samplecov');
        pe_blockdiag =  xlsread(intfname,'MAIN','pe_blockdiag');
        pe_blockdiagexc =  xlsread(intfname,'MAIN','pe_blockdiagexc');
        [junk pe_country_exc] =  xlsread(intfname,'MAIN','pe_exccountry'); %#ok
        
        if pe_samplecov == 1
            % use original estimated varcov matrix
            pe_meth = 1;
        elseif pe_blockdiag == 1
            % transform the original varcov matrix to a block diagonal one
            pe_meth = 2;
        elseif pe_blockdiagexc == 1
            % transform the original varcov matrix to a block diagonal one,
            % while leaving unrestricted the cross-country correlations for a
            % given exception country
            pe_meth = 3;
        end
        
        % Settings for shrinkage
        %**********************************************************************
        % note: the varcov matrix is used only for data generation in this
        % section
        
        use_shrinkedvcv_dg =  xlsread(intfname,'MAIN','use_shrinkedvcv_dg');
        
        % Recover info about the type of bootstrap
        [junk Bootstrap_approach] = xlsread(intfname,'MAIN','Bootstrap_approach');
        
        if strcmp('inverse',Bootstrap_approach) == 1
            shuffleflag=0;
        elseif strcmp('shuffle',Bootstrap_approach) == 1
            shuffleflag=1;
        end
        
        
        % Transform the varcov matrix for bootstrap:
        pe_varcov_tx = transform_varcov(pe_meth,pe_country_exc,Sigma_zeta,cnames_s_y);
        
        
        
        % Shrinkage of variance covariance matrix for bootstrap
        %**************************************************************
        arbitrarylambda = []; do_alambda = 0;
        if pauseflag == 0 % no pauses
            % check if the user has inputed an arbitrary value for the
            % shrinkage parameter lambda
            arbitrarylambda = xlsread(intfname,'MAIN','arbitrarylambda');
            if not(isempty(arbitrarylambda))
                do_alambda = 1;
            end
        end
        
        
        if shuffleflag == 1 % shuffle-option of bootstrap
            % no shrinkage
            pe_varcov_dg = pe_varcov_tx;
            lambda_star = [];
        elseif shuffleflag == 0 % inverse-option of bootstrap, % need to do shrinkage
            
            % check if the pe_varcov_tx matrix is positive definite
            posflag=1;
            try
                chol(pe_varcov_tx);
            catch %#ok
                posflag=0;
            end
            firstpause = 1;
            
            if posflag == 1
                if use_shrinkedvcv_dg == 0
                    % don't do shrinkage
                    pe_varcov_dg=pe_varcov_tx;
                elseif use_shrinkedvcv_dg == 1
                    if do_alambda == 1
                        % shrinkage
                        lambda_star = arbitrarylambda;
                        pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                    else
                        [pe_varcov lambda_star]=ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                        
                        % print the resulting shrinkage parameter in the
                        % interface file
                        xlswrite(intfname,lambda_star,'MAIN','EH33');
                        
                    end
                end
            elseif posflag == 0
                k=1;
                while k<2
                    
                    disp('>>> Warning: The covariance matrix is not positive definite.')
                    disp(' ');
                    if firstpause == 1 || (firstpause == 0 && isempty(arbitrarylambda))
                        msg = sprintf('>>> Pause and go to %s: Ensure that perform shrinkage on the correlation matrix',intfname);
                        disp(msg);
                        disp('        for generating the bootstrap data is set to 1 (if it is not already),');
                        disp('        then press enter.');
                    elseif firstpause == 0 && not(isempty(arbitrarylambda))
                        msg = sprintf('>>> Pause and go to %s: Increase the value of the shrinkage parameter',intfname);
                        disp(msg);
                        disp('         without making any other changes, then press enter.');
                    end
                    disp(' ');
                    winopen(intfname);
                    pause
                    msg = sprintf('>>> Make sure you have saved and closed the %s file. If so, press enter again.',intfname);
                    disp(msg);
                    disp('The program is now running (do not press any key)  ');
                    disp(' ');
                    
                    arbitrarylambda =  xlsread(intfname,'MAIN','arbitrarylambda');
                    use_shrinkedvcv_dg =  xlsread(intfname,'MAIN','use_shrinkedvcv_dg');
                    firstpause = 0;
                    
                    
                    if use_shrinkedvcv_dg == 1
                        if do_alambda == 1
                            lambda_star = arbitrarylambda;
                            pe_varcov_dg =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),lambda_star);
                        else
                            [pe_varcov_dg lambda_star] =ShrinkageCorrLstar(pe_varcov_tx,cols(zeta),[]);
                            k = 2;
                        end
                    elseif use_shrinkedvcv_dg == 0
                        pe_varcov_dg = pe_varcov_tx;
                    end
                    
                    posflag=1;
                    try
                        chol(pe_varcov);
                    catch %#ok
                        posflag=0;
                    end
                    
                    if posflag == 1
                        k = 2;
                        disp('- OK: The covariance matrix is positive definite');
                    end
                end
            end
        end
        
        
        
        [kpsup_90cv kpmnsq_90cv ny_90cv rny_90cv qlr_90cv mw_90cv apw_90cv rqlr_90cv rmw_90cv rapw_90cv...
            kpsup_95cv kpmnsq_95cv ny_95cv rny_95cv qlr_95cv mw_95cv apw_95cv rqlr_95cv rmw_95cv rapw_95cv ...
            kpsup_99cv kpmnsq_99cv ny_99cv rny_99cv qlr_99cv mw_99cv apw_99cv rqlr_99cv rmw_99cv rapw_99cv] = bootstrap_GVAR_ss(ss_B,...
            ccut,kpsup,kpmnsq,ny,rny,qlr,mw,apw,rqlr,rmw,rapw,zeta, ...
            pe_varcov_dg,y,ynames,cnum,cnames,...
            delta_0,delta_1,C,H0,dvtype,ntypes,wmatrices,maxlag,mlag,varxlag,estcase,rank,...
            dv, vnames, vnames_long, fvnames, gvnames, gvnames_long, vnum, gvnum, gxidx, gnvnum,gnvnames, ...
            dvflag, fvflag, gvflag,ydate,ylabelseq_estimation,yearseq_estimation,nyears_estimation,nobs,overid_flag,beta_r,shuffleflag);
        
        if printout == 1
            disp('- Adding to output.xlsx: Structural stability test, bootstrapped critical values');
            print_sstests_cvals(kpsup_90cv,kpmnsq_90cv,ny_90cv,rny_90cv,qlr_90cv,mw_90cv,apw_90cv,rqlr_90cv,rmw_90cv,rapw_90cv,...
                kpsup_95cv,kpmnsq_95cv,ny_95cv,rny_95cv,qlr_95cv,mw_95cv,apw_95cv,rqlr_95cv,rmw_95cv,rapw_95cv,...
                kpsup_99cv,kpmnsq_99cv,ny_99cv,rny_99cv,qlr_99cv,mw_99cv,apw_99cv,rqlr_99cv,rmw_99cv,rapw_99cv,...
                cnum,cnames,cnames_long,outdir);
        end
    end
end


disp('--- End of the program ---'); 

