function rlevelres = aggregate_results(rweights,vnames,gvnames,rnames,...
    regions,xnames,cnames_s_x,clevelres)

%**************************************************************************
% PURPOSE: Aggregate GIRF or GFEVD country-level figures into 
%          regional-level ones 
%--------------------------------------------------------------------------
% INPUT:
% - rweights: struct, for each variable it contains weights of each country
% wrt to the region to which it belongs 
% - vnames: cell, list of domestic variables
% - gvnames: cell,list of global variables 
% - rnames: cell, list of regions for aggregation of results 
% - regions: struct, contains lists of countries of each region 
% - xnames: cell, ordered list of short names of each endogenous variable 
% in the GVAR (has dimension K times 1)
% - cnames_s_x: cell, ordered list of country-labels (short names) of each 
% endogenous variable in the GVAR (has dimension K times 1) 
% - clevelres: double, matrix with dimension K x (Forecast Horizon + 1)
% containing country-level GIRFs or GFEVDs
%--------------------------------------------------------------------------
% OUTPUT: 
% - rlevelres: struct, contains region-level GIRFs or GFEVDS for each
% variable 
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

allvnames = [vnames gvnames];

for k=1:length(allvnames)
    for r=1:length(rnames)
        for c=1:length(regions.(rnames{r}))
            for j=1:rows(clevelres)
                if (strcmp(xnames{j},allvnames{k}) && strcmp(cnames_s_x{j},regions.(rnames{r}){c}))
                    clevelres_str.(allvnames{k}).(rnames{r}).(regions.(rnames{r}){c}) = clevelres(j,:);
                end
            end
        end
    end
end

for k=1:length(allvnames)
    for r=1:length(rnames)
        rweightsblock = [];
        gblock = [];

        for c = 1:length(regions.(rnames{r}))
            if isfield(rweights,allvnames{k})
                if isfield(rweights.(allvnames{k}),rnames{r}) && isfield(rweights.(allvnames{k}).(rnames{r}),regions.(rnames{r})(c))
                    rweightsblock = [rweightsblock rweights.(allvnames{k}).(rnames{r}).(regions.(rnames{r}){c})]; %#ok
                end
            end
            if isfield(clevelres_str,allvnames{k})
                if isfield(clevelres_str.(allvnames{k}),rnames{r}) && isfield(rweights.(allvnames{k}).(rnames{r}),regions.(rnames{r})(c))
                    gblock = [gblock; clevelres_str.(allvnames{k}).(rnames{r}).(regions.(rnames{r}){c})]; %#ok
                end
            end
        end
        if isfield(clevelres_str,allvnames{k})
            if isfield(clevelres_str.(allvnames{k}),rnames{r})
                rlevelres.(allvnames{k}).(rnames{r}) = rweightsblock*gblock;
            end
        end
    end
end
            