
function [endog endoglist endoglist_long exog exoglist exognx exognxlist] = create_countrymodels(cnames,cnum, ...
    vnames,vnames_long,fvnames,vnum,gvnames,gvnames_long,gvnum,gnvnum,gnvnames,dv,fv,gv,dvflag,fvflag,gvflag)

%**************************************************************************
% PURPOSE: create individual models (e.g. country models)
%--------------------------------------------------------------------------
% INPUT: 
% - cnames: list of countries' names (short names)
% - cnum: number of countries
% - vnames: contains the list of domestic variables (short names)
% - vnames_long: contains the list of domestic variables (long names)
% - fvnames: contains the list of foreign-specific variables (short names)
% - vnum: number of domestic variables
% - gvnames: contains the list of global variables (short names)
% - gvnames_long: contains the list of global variables (long names)
% - gvnum: number of global variables
% - gnvnum: number of global variables which are endogenous in the GVAR
% - gnvnames: contains the list of global variables which are endogenous 
%   in the GVAR(short names)
% - dv: structure, contains data of domestic variables for each country
% - fv: structure, contains data of foreign-specific variables for each 
%   country
% - gv: structure, contains data of global variables 
% - dvflag: cnum x vnum matrix of dummies obtained from gvar.xls,
%   indicating user's choice of domestic variables for each country
% - fvflag: cnum x fvnum matrix of dummies obtained from gvar.xls,
%   indicating user's choice of foreign-specific variables for each country
% - gvflag: cnum x gvnum matrix of dummies obtained from gvar.xls,
%   indicating user's choice of global variables for each country
%--------------------------------------------------------------------------
% OUTPUT: 
% - endog: structure, contains data of endogenous variables of each country
% - endoglist: structure, contains list of endogenous variables of each
%   country (short names)
% - endoglist_long: structure, contains list of endogenous variables of each
%   country (long names)
% - exog: structure, contains data of weakly exogenous variables of each 
%   country
% - exoglist: structure, contains list of weakly exogenous variables of 
%   each country (short names)
% - exognx: structure, contains data of weakly exogenous variables of each 
%   country, excluding global variables which are exogenous to the GVAR
% - exognxlist: structure, contains list of weakly exogenous variables of 
%   each country, excl global variables which are exogenous to the GVAR (short names)
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


for n=1:cnum
    endog.(cnames{n}) = [];   % creating the block of endogenous variables
    endoglist.(cnames{n}) = {}; % creating the list of endogenous var names
    endoglist_long.(cnames{n}) = {}; % creating the list of endogenous var full names
    
    exog.(cnames{n}) = [];  % creating the block of (weakly) exogenous variables
    exoglist.(cnames{n}) = {}; % creating the list of exogenous var names
    
    exognx.(cnames{n}) = []; % creating the block of (weakly) exogenous variables, excluding
    % the global variables which are exogenous to the GVAR
    exognxlist.(cnames{n}) = {}; % creating the corresponding list of var names
    
    for i=1:vnum
        if dvflag(n,i) == 1
            endog.(cnames{n}) = [endog.(cnames{n}) dv.(vnames{i}).(cnames{n})];
            endoglist.(cnames{n}) = [endoglist.(cnames{n}) vnames{i}];
            endoglist_long.(cnames{n}) = [endoglist_long.(cnames{n}) vnames_long{i}];
        end
          
        if fvflag(n,i) == 1
            exog.(cnames{n}) = [exog.(cnames{n}) fv.(fvnames{i}).(cnames{n})];
            exoglist.(cnames{n}) = [exoglist.(cnames{n}) fvnames{i}];   
            
            exognx.(cnames{n}) = exog.(cnames{n});
            exognxlist.(cnames{n}) = exoglist.(cnames{n});
        end
    end
    

    if not(gvnum==0)
        % adding the global variables
        for i=1:gvnum

            if gvflag(n,i) == 0 % global variable does not enter in the model

            elseif gvflag(n,i) == 1  % global variable is weakly exogenous
                exog.(cnames{n}) = [exog.(cnames{n}) gv.(gvnames{i})];
                exoglist.(cnames{n}) = [exoglist.(cnames{n}) gvnames{i}];
                
                globvarname = gvnames(i);
                if not(gnvnum==0) % if there are global variables which are endogenous
                    for j=1:gnvnum
                        if strcmp(gnvnames(j),globvarname)
                            exognx.(cnames{n}) = [exognx.(cnames{n}) gv.(gvnames{i})];
                            exognxlist.(cnames{n}) = [exognxlist.(cnames{n}) gvnames{i}];
                        end
                    end
                end
            elseif gvflag(n,i) == 2 % global variable is endogenous
                endog.(cnames{n}) = [endog.(cnames{n}) gv.(gvnames{i})];
                endoglist.(cnames{n}) = [endoglist.(cnames{n}) gvnames{i}];
                endoglist_long.(cnames{n}) = [endoglist_long.(cnames{n}) gvnames_long{i}];
            end
        end
    end
end
       


