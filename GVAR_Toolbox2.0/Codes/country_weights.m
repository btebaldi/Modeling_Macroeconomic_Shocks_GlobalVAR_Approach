
function cweights = country_weights(aggrwgts,cnames,vnum,vnames,...
    gvnum,gvnames,endoglist,wtypes)
                                    
%**************************************************************************
% PURPOSE: calculates country weights for each domestic variable, which are
% used for computing global shocks and for computing regional weights
%--------------------------------------------------------------------------
% INPUT:
% - aggrwgts: struct, it contains country-level data for computing 
% country weights (e.g. PPP-GDP figures) 
% - cnum: number of country models
% - cnames: cell, list of countries' names (short names)
% - vnum: number of domestic variables
% - vnames: cell, list of domestic variables (short names) 
% - gvnum: number of global variables
% - gvnames: cell, list of global variables (short names)
% - endoglist: struct, it contains lists of endogenous variables for each
% country
%--------------------------------------------------------------------------
% OUTPUT:
% - cweights: struct, for each variable it contains countries' weights 
%--------------------------------------------------------------------------
% NOTES:  reweighting of weights if for a given country a given variable is 
%         missing
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


cnum = length(cnames);

    % for creating the feedback variables
    ntypes = max(wtypes);
    all_aggrwgts = zeros(cnum,ntypes);
    for n=1:cnum
         rsaggrwgts=aggrwgts.(cnames{n});
        all_aggrwgts(n,:) = rsaggrwgts(1:ntypes);
    end


% make weights that sum to one for each domestic variable

flag = zeros(vnum,cnum);
for j=1:vnum
    if wtypes == 1
        tmp = all_aggrwgts;
    else
        weighttype = wtypes(j);
        if weighttype == 0 % the variable j is not going to be used as feedback variable
            weighttype = 1; % select weight type = 1, anyway the variable is not used
        end
        tmp = all_aggrwgts(:,weighttype);
    end
    
    
    for n=1:cnum
        
        
        for k=1:length(endoglist.(cnames{n}))
            if strcmp(vnames(j),endoglist.(cnames{n})(k))
                flag(j,n) = 1;
            end
        end
        
        if flag(j,n) == 0
            tmp(n) = 0;
        end
    end
    
    for n=1:cnum
        cweights_tmp = tmp(n)./sum(tmp);
        if not(flag(j,n) == 0)
            cweights.(vnames{j}).(cnames{n}) = cweights_tmp;
        end
    end
end

if wtypes == 1
    if not(gvnum==0)
        gflag = zeros(gvnum,cnum);
        for g=1:gvnum
            for n=1:cnum
                for j=1:length(endoglist.(cnames{n}))
                    if strcmp(gvnames(g),endoglist.(cnames{n})(j))
                        gflag(g,n) = 1;
                    end
                end
                if gflag(g,n) == 0
                    tmp(n) = 0;
                else
                    tmp(n) = 1;
                end
            end
            
            for n=1:cnum
                cweights_tmp = tmp(n);
                if not(gflag(g,n) == 0)
                    cweights.(gvnames{g}).(cnames{n}) = cweights_tmp;
                end
            end
        end
    end
end


