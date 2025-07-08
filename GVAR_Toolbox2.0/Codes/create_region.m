function [dv aggrwgts] =  create_region(name,list,nobs,vnum,vnames,dv,aggrwgts)

%**************************************************************************
% PURPOSE: aggregate country-specific data and weights in order to form
% regional data
%--------------------------------------------------------------------------
% INPUT:
% - name: cell object, contains the name of the region that will be 
% constructed
% - list: cell object, contains the list of countries belonging to region 
% 'name'
% - nobs: number of observations
% - vnum: number of domestic variables
% - vnames: cell object, contains the list of domestic variables
% - dv: structure, contains data of domestic variables for each country
% - aggrwgts: structure, contains weights (note they not yet add up to one) 
%             for aggregation of all the countries
%--------------------------------------------------------------------------
% OUTPUT:
% - dv: structure, updated as now contains data of the brand new region and
% does not contain data of corresponding individual countries anymore
% - aggrwgts: structure, updated as now contains weight of the brand new
% region in place of individual weights of individual countries belonging
% to the region.
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************


aggrwgtsblock=[];
for j=1:length(list)
    
    aggrwgtsblock = [aggrwgtsblock aggrwgts.(list{j})]; %#ok
    
    aggrwgts = rmfield(aggrwgts, list{j});
    
end

aggrwgts.(name{1}) = sum(aggrwgtsblock);
aggrwgts = orderfields(aggrwgts);

for i=1:vnum
    
    regweights = aggrwgtsblock/sum(aggrwgtsblock);
    
    
    varblock=[];
    for j=1:length(list)
        
        
        
        if isfield(dv.(vnames{i}),list{j})
        varblock = [varblock dv.(vnames{i}).(list{j})]; %#ok
           
        % % delete region countries data        
        dv.(vnames{i}) = rmfield(dv.(vnames{i}), list{j});
        else
        varblock = [varblock zeros(nobs,1)]; %#ok
        regweights(j) = 0;
        end
        
        
        
    end
    
    % resize weights
    regweights = regweights./sum(regweights);
    
    dv.(vnames{i}).(name{1}) = varblock*regweights';
    
    
    % check whether one region has not a particular variable, if so, remove
    % the corresponding field
    if isnan(dv.(vnames{i}).(name{1})(1)) % missing that variable
        dv.(vnames{i}) = rmfield(dv.(vnames{i}),name{1});
    end
    
    % reorder structure
    dv.(vnames{i}) = orderfields(dv.(vnames{i}));
end
      
