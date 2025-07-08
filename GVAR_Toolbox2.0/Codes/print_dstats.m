function print_dstats(vnum,vnames,vnames_long,cnum,cnames,cnames_long,dv,...
    fv,gv,gvnum,gvnames,gvnames_long,outdir)

%**************************************************************************
% PURPOSE: printing descriptive statistics of domestic, foreign-specific
% and global variables
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     

% Domestic variables
%**************************************************************************
tab = ['Descriptive Statistics of Domestic Variables' num2cell(NaN(1,9)); ' ' num2cell(NaN(1,9))];
for j =1:vnum
    tab2 = {vnames_long{j} 'Mean' 'Median' 'Maximum' 'Minimum' 'Std. dev.' 'Skewness' 'Kurtosis' 'Jarque-Bera' 'Probability'};
    tab = [tab; tab2; '*******************' num2cell(NaN(1,9))]; %#ok
    
    for n=1:cnum
        if isfield(dv.(vnames{j}),cnames{n})
            [mn md mx minv st sk ku] = dstats(dv.(vnames{j}).(cnames{n}));
            [W chi2_pval] = jarquebera(dv.(vnames{j}).(cnames{n}));
            tab = [tab; cnames_long{n} num2cell([mn md mx minv st sk ku W chi2_pval])]; %#ok
        else
            tab = [tab; cnames_long{n} num2cell(NaN(1,9))]; %#ok
        end
    end
    tab = [tab; ' ' num2cell(NaN(1,9))]; %#ok
end
xlswrite([outdir 'output.xlsx'],tab,'dstats_domestic');


% Foreign-specific variables
%**************************************************************************
tab = ['Descriptive Statistics of Foreign-Specific Variables' num2cell(NaN(1,9)); ' ' num2cell(NaN(1,9))];
for j =1:vnum
    tab2 = {vnames_long{j} 'Mean' 'Median' 'Maximum' 'Minimum' 'Std. dev.' 'Skewness' 'Kurtosis' 'Jarque-Bera' 'Probability'};
    tab = [tab; tab2; '*******************' num2cell(NaN(1,9))]; %#ok
    
    for n=1:cnum
        if isfield(fv.(vnames{j}),cnames{n})
            [mn md mx minv st sk ku] = dstats(fv.(vnames{j}).(cnames{n}));
            [W chi2_pval] = jarquebera(fv.(vnames{j}).(cnames{n}));
            tab = [tab; cnames_long{n} num2cell([mn md mx minv st sk ku W chi2_pval])]; %#ok
        else
            tab = [tab; cnames_long{n} num2cell(NaN(1,9))]; %#ok
        end
    end
    tab = [tab; ' ' num2cell(NaN(1,9))]; %#ok
end
xlswrite([outdir 'output.xlsx'],tab,'dstats_foreign');



if not(gvnum==0)
    % Global variables
    %**********************************************************************
    tab = ['Descriptive Statistics of Global Variables' num2cell(NaN(1,9)); ' ' num2cell(NaN(1,9))];
    tab2 = {'Statistics' 'Mean' 'Median' 'Maximum' 'Minimum' 'Std. dev.' 'Skewness' 'Kurtosis' 'Jarque-Bera' 'Probability'};
    tab = [tab; tab2; '*******************' num2cell(NaN(1,9))];
    
    for j=1:gvnum
        [mn md mx minv st sk ku] = dstats(gv.(gvnames{j}));
        [W chi2_pval] = jarquebera(gv.(gvnames{j}));
        tab = [tab; gvnames_long{j} num2cell([mn md mx minv st sk ku W chi2_pval])]; %#ok
    end
    tab = [tab; ' ' num2cell(NaN(1,9))];
    xlswrite([outdir 'output.xlsx'],tab,'dstats_global');
end


    

