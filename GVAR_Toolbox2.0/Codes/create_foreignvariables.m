function fv = create_foreignvariables(cnum,cnames,vnum,vnames,dv,dvtype,wmatrices)

%**************************************************************************
% PURPOSE: create foreign-specific variables
%--------------------------------------------------------------------------
% INPUT: 
% - cnum: number of countries
% - cnames: list of countries' names (short names)
% - vnum: number of domestic variables
% - vnames: cell object, contains the list of domestic variables
% - dv: structure, contains data of domestic variables for each country
% - dvtype: structure, contains info about type of each variable 
% - wmatrices: structure containing n weight matrices (n is number of 
%              variable types) 
%--------------------------------------------------------------------------
% OUTPUT: 
% - fv: structure, contains data of foreign-specific variables for each
% country
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************



for i = 1:vnum
    
    weight_nature = dvtype.(vnames{i});
    wmattype = sprintf('wmat%d',weight_nature);
    weightsmatrix = wmatrices.(wmattype);

    vlist = fieldnames(dv.(vnames{i}));
    vn = length(vlist);
    
    vblock = [];
    for n=1:vn
        vblock = [vblock dv.(vnames{i}).(vlist{n})]; %#ok
    end
    
    vind = [];
    for n=1:cnum
        for j=1:vn
            if (strcmp(cnames{n},vlist{j}))
                vind = [vind; n]; %#ok
            end
        end
    end
    
    % resize weight matrix in order to skip countries without that variable
    wmatx.(vnames{i}) = weightsmatrix(vind,:);
    
    % reweight weightsmatrix in order that all columns sum to one
    for j=1:cols(wmatx.(vnames{i}))
        wmatx.(vnames{i})(:,j) = wmatx.(vnames{i})(:,j)./sum(wmatx.(vnames{i})(:,j));
    end
    
    % create foreign variables
    fvblock = vblock*wmatx.(vnames{i}); % block of foreign variables
    
    % create structure
    for n=1:cnum
        fv.(vnames{i}).(cnames{n}) = fvblock(:,n);
    end
    
end
