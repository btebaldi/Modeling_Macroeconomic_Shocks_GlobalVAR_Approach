function print_rweights(vnames,gvnames,rnames,regions,cweights,rweights,outdir)

%**************************************************************************
% PURPOSE: printing regional weights for aggregation
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     

title = {'Regional Weights'};
hlabel = ['Region' 'Country' vnames gvnames];
vlabel = [];
for r=1:length(rnames)
    for j=1:length(regions.(rnames{r}))
        vlabel = [vlabel; rnames(r) regions.(rnames{r})(j)]; %#ok
    end
end
table = [];
for j=1:length(vnames)
    wgt_tmp=[];
    for r=1:length(rnames)
        for n=1:length(regions.(rnames{r}))
            if isfield(cweights.(vnames{j}),regions.(rnames{r}){n});
                wgt_tmp = [wgt_tmp; rweights.(vnames{j}).(rnames{r}).(regions.(rnames{r}){n})]; %#ok
            else
                wgt_tmp = [wgt_tmp; NaN]; %#ok
            end
        end
    end
    table = [table wgt_tmp]; %#ok

end

for g=1:length(gvnames)
    wgt_tmp = [];
    for r=1:length(rnames)
        for n=1:length(regions.(rnames{r}))
            if isfield(cweights,gvnames{g})
                if isfield(cweights.(gvnames{g}),regions.(rnames{r}){n});
                    wgt_tmp = [wgt_tmp; rweights.(gvnames{g}).(rnames{r}).(regions.(rnames{r}){n})]; %#ok
                else
                    wgt_tmp = [wgt_tmp; NaN]; %#ok
                end
            end
        end
    end
    table = [table wgt_tmp]; %#ok
end

tab = [title num2cell(NaN(1,1+cols(table)))];
tab = [tab; num2cell(NaN(1,2+cols(table)))];
tab = [tab; hlabel; vlabel num2cell(table)];

xlswrite([outdir 'output.xlsx'],tab,'rwgts');

        
