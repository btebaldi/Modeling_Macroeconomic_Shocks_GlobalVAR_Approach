function print_VARXord(cnames,cnames_long,cnum,varxlag,outdir)


%**************************************************************************
% PURPOSE: printing lag orders of VARX* models
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
tab = ['VARX* Order of Individual Models (p: lag order of domestic variables, q: lag order of foreign variables)' num2cell(NaN(1,2))];
tab = [tab; num2cell(NaN(1,3)); {' ' 'p' 'q'}];

out = [];
for n=1:cnum
    out = [out; cnames_long(n) num2cell(varxlag.(cnames{n}))]; %#ok
end

tab = [tab; out];

xlswrite([outdir 'output.xlsx'],tab,'VARXorder');


