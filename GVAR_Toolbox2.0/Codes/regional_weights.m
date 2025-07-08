
function rweights = regional_weights(vnames,gvnames,rnames,regions,cweights)

%**************************************************************************
% PURPOSE: calculates regional weights for each domestic variable, for
% regional aggregation of shocks and results
%--------------------------------------------------------------------------
% INPUT:
% - vnames: cell, list of domestic variables (short names) 
% - gvnames: cell, list of global variables (short names)
% - rnames: list of regions for aggregation (short names) 
% - regions: struct, contains lists of countries of each region 
% - cweights: struct, for each variable it contains countries' weights 
%--------------------------------------------------------------------------
% OUTPUT:
% - rweights: struct, for each variable it contains regional weights 
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
allvnames = [vnames gvnames];
allvnum = length(allvnames);

rnum = length(rnames);

for i=1:allvnum
    for r=1:rnum
        sumcweights = 0;
        for j=1:length(regions.(rnames{r}))
            if isfield(cweights,allvnames{i})
                if isfield(cweights.(allvnames{i}),regions.(rnames{r}){j});
                    sumcweights = sumcweights + cweights.(allvnames{i}).(regions.(rnames{r}){j});
                end
            end
        end
        for j=1:length(regions.(rnames{r}))
            if isfield(cweights,allvnames{i})
                if isfield(cweights.(allvnames{i}),regions.(rnames{r}){j});
                    rweights.(allvnames{i}).(rnames{r}).(regions.(rnames{r}){j}) = cweights.(allvnames{i}).(regions.(rnames{r}){j})/sumcweights;
                end
            end
        end
    end
end
