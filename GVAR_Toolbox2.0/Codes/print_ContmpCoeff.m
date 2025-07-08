function print_ContmpCoeff(cnum,cnames,cnames_long,ctplabel,ctpcoeffs,ctpstd,...
         ctpstd_hcw,ctpstd_nwc,ctptvals,ctptvals_hcw,ctptvals_nwc,vnames,outdir)

%**************************************************************************
% PURPOSE: printing contemporaneous effects of the foreign variables on
% their corresponding domestic counterparts
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
  
tab = ['Contemporaneous Effects of Foreign Variables on Domestic Counterparts' num2cell(NaN(1,1+length(vnames)))];
tab = [tab; num2cell(NaN(1,2+length(vnames)))];
tab = [tab; num2cell(NaN(1,2)) vnames];

vlabel1 = [];
vlabel2_tmp = {'Coefficient'; 'Standard error'; 't-ratio'; ...
    'White''s adjusted SE'; 't-ratio_White'; 'Newey-West''s adjusted SE'; 't-ratio_NeweyWest'};
vlabel2 = [];
table = [];

k=1;
for n=1:cnum
    if isfield(ctplabel,cnames(n)) % skip country-models which have no foreign-specific variables
        table = [table; NaN(7,length(vnames))];
        
        clabel = cnames_long(n); %#ok
        for z=1:7
            vlabel1 = [ vlabel1; cnames_long(n)]; %#ok
        end
        
        vlabel2 = [vlabel2; vlabel2_tmp]; %#ok
        for i=1:length(ctplabel.(cnames{n}))
            for j=1:length(vnames)
                if strcmp(ctplabel.(cnames{n})(i),vnames(j))
                    table(k,j) = ctpcoeffs.(cnames{n})(i);
                    table(k+1,j) = ctpstd.(cnames{n})(i);
                    table(k+2,j) = ctptvals.(cnames{n})(i);
                    table(k+3,j) = ctpstd_hcw.(cnames{n})(i);
                    table(k+4,j) = ctptvals_hcw.(cnames{n})(i);
                    table(k+5,j) = ctpstd_nwc.(cnames{n})(i);
                    table(k+6,j) = ctptvals_nwc.(cnames{n})(i);
                else
                    continue
                end
            end
        end
        k=k+7;
    end
end

tab = [tab; vlabel1 vlabel2 num2cell(table)];

xlswrite([outdir 'output.xlsx'],tab,'contemp_coeff');

