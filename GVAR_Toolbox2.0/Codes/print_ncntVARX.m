function print_ncntVARX(cnum,cnames,cnames_long,rank,outdir)

%**************************************************************************
% PURPOSE: printing ranks of each individual VECMX* model
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     


tab = ['# Cointegrating Relationships for the Individual VARX* Models' num2cell(NaN(1,1))];
tab = [tab; num2cell(NaN(1,2))];
tab = [tab; {'Country' '# Cointegrating relations'}];

rank_out = [];

for n=1:cnum
    rank_out = [rank_out; rank.(cnames{n})]; %#ok
end
    
tab = [tab; cnames_long num2cell(rank_out)];


xlswrite([outdir 'output.xlsx'],tab,'ncntVARX');
