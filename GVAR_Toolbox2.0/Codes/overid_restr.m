function [beta_r nunrestrpar] = overid_restr(cnum,cnames,cnames_long,endoglist,...
         exoglist,gvnum,gvnames,estcase,outdir,rank)

%**************************************************************************
% PURPOSE: It allows the user to input overidentifying restrictions, then
% yields beta_r (the imputed restrictions) and nunrestrpar, the number of
% unrestricted coefficients in all the cointegrating relations of each
% country model
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
      
 


% creating worksheet in gvar.xls for inputing restrictions
%*********************************************************

fcells = 100;

tab = ['Insert Overidentifying Restrictions on the Cointegrating Vectors' num2cell(NaN(1,fcells-1))];
tab = [tab; num2cell(NaN(1,fcells))];

pointer = 3;
pointerlist = zeros(cnum,1);

for n=1:cnum
    vlabel = cnames_long(n); 
    hlabel = {};
    
    endlist = endoglist.(cnames{n});
    exolist_t = exoglist.(cnames{n});
    
    exolist = {};
    for j=1:length(exolist_t)
        flag = 0;
        if not(gvnum==0)
            for g=1:gvnum
                if strcmp(exolist_t(j),gvnames(g))
                    flag = 1;
                end
            end
        end
        
        if flag == 0
            tmp = sprintf('%ss',exolist_t{j});
        else
            tmp = sprintf('%s',exolist_t{j});
        end
        tmp = {tmp};
        exolist = [exolist tmp]; %#ok
    end
    
    if estcase.(cnames{n}) == 2
        hlabel = ['Intercept' endlist exolist];
    elseif estcase.(cnames{n}) == 3
        hlabel = [endlist exolist];
    elseif estcase.(cnames{n}) == 4
        hlabel = ['Trend' endlist exolist];
    end
    
    tab = [tab; vlabel hlabel num2cell(NaN(1,fcells-length(vlabel)-length(hlabel)))]; %#ok
    
    
    
    cbeta = rank.(cnames{n});
    
    
    for j=1:cbeta
        tmp = sprintf('CV%d',j);
        tmp = {tmp};
        tab = [tab; tmp num2cell(NaN(1,fcells-1))]; %#ok
    end
    tab = [tab; '# unrestricted:' num2cell(NaN(1,fcells-1))]; %#ok
    
    
    if n>1
        cbetap = rank.(cnames{n-1});
        
        pointer = pointer + (cbetap+1) +1;
    end
    
    pointerlist(n) = pointer +1;
end

xlswrite([outdir 'overid_restr.xls'],tab,'overid_restr');


% pause and let user input overidentifying restrictions
%**********************************************************************
winopen([outdir 'overid_restr.xls']);
disp(' ');
disp('>>> Pause and go to overid_restr.xls: Insert overidentifying restrictions on the');
disp('    cointegrating vector(s) of the individual VECMX* models, then press enter.');
disp(' ');
disp('=========================================================================================');
disp('Note that:');
disp('1. Overidentifying restrictions should be imposed on the coefficients of ALL cointegrating');
disp('   relations of a particular country, though not necessarily on all countries. You only');
disp('   need to fill the cells corresponding to the cointegrating vector(s) of those countries');
disp('   (or country) for which you wish to impose restrictions.');
disp('2. Restrictions should be imposed on ALL elements of any cointegrating vector. If you');
disp('   wish to allow any of the elements to be unrestricted, these have to be estimated outside');
disp('   the program and subsequently imposed here. The total number of these elements, across all'); 
disp('   cointegrating vectors for a particular country, should be imposed in the cell adjacent');
disp('   to # unrestricted, otherwise this cell should be left empty. See the user guide for');   
disp('   further information.');
disp('4. The cointegrating vectors associated with empty cells will be estimated by default');
disp('   under exact identification, using the identity matrix normalisation scheme.')
disp('5. All restrictions inserted, will be imposed simultaneously during estimation.');   
disp(' ');
disp(' USING THE FULL DEMO INTERFACE FILE');
disp('If you are using the full demo interface file, and would like to replicate the results');
disp('in the Output Full Demo folder downloaded with the program:                           ');
disp('1. For Canada:');
disp('CV1:  Set Dp=-1, r=1, and zero for the rest of the variables');
disp('CV2:  Set lr=1, lrs=-1, and zero for the rest of the variables');
disp('CV3:  Set r=-1, lr=1, and zero for the rest of the variables');
disp('2. For Euro:');
disp('CV1:  Set Dp=-1, r=1, and zero for the rest of the variables');  
disp('CV2:  Set r=-1, lr=1, and zero for the rest of the variables');
disp('=========================================================================================');
disp(' ');
pause
disp('>>> Make sure you have saved and closed the overid_restr.xls file. If so, press enter again.');
endpause = input(' ','s'); %#ok
disp('- The program is now running (do not press any key)');


% read inputed overidentifying restrictions
%********************************************

restrlist = zeros(cnum,1);
restridx = [];
for n=1:cnum
    
    location = sprintf('B%d',pointerlist(n));
    
    restr_t = xlsread([outdir 'overid_restr.xls'],'overid_restr',location);
    
    if not(isnan(restr_t))
        restrlist(n) = 1;
        restridx = [restridx; n]; %#ok
    end
    
end

nrestr = sum(restrlist);
        
[allrestr_t junk] = xlsread([outdir 'overid_restr.xls'],'overid_restr'); %#ok

allrestr = [];
nunrestr = zeros(nrestr,1); k = 1;
for i=1:rows(allrestr_t)
    if not(isnan(allrestr_t(i,1:2)))
        allrestr = [allrestr; allrestr_t(i,:)]; %#ok
    end
    
    if not(isnan(allrestr_t(i,1))) && isnan(allrestr_t(i,2))
        nunrestr(k) = allrestr_t(i,1);
        k = k+1;
    elseif (isnan(allrestr_t(i,1)) && isnan(allrestr_t(i,2))) && not(isnan(allrestr_t(i-1,1))) && not(isnan(allrestr_t(i-1,2)))
        nunrestr(k) = 0;
        k = k+1;
    end
    
end

if nrestr>0
    k=1;
    for n=1:nrestr
        cbeta = rank.(cnames{restridx(n)});
        beta_r_tmp = allrestr(k:k+cbeta-1,:);
        beta_rt = [];
        for j=1:cols(beta_r_tmp)
            if not(isnan(beta_r_tmp(1,j)))
                beta_rt = [beta_rt beta_r_tmp(:,j)]; %#ok
            end
        end

        beta_r.(cnames{restridx(n)}) = beta_rt';
        nunrestrpar.(cnames{restridx(n)}) = nunrestr(n);
        
        k=k+cbeta;
    end
else
    beta_r = [];
    nunrestrpar = [];
end

