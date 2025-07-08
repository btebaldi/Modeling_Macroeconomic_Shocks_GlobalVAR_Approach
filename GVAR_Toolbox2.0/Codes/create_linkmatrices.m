
function [z x xnames W] = create_linkmatrices(cnames,cnum,vnum,vnames,...
    fvnum,fvnames,gvnum,gvnames,endog,endoglist,exog,exoglist,dvtype,wmatrices)

%**************************************************************************
% PURPOSE: builds the countries' series z_it, the global vector x_t, and 
%          then computes the link matrices W_i of each country i
%--------------------------------------------------------------------------
% INPUT:
% - cnames: list of countries' names (short names)
% - cnum: number of countries
% - vnum: number of domestic variables
% - vnames: cell, contains the list of domestic variables (short names)
% - fvnum: number of foreign-specific variables
% - fvnames: cell, contains the list of domestic variables (short names)
% - gvnum: number of global variables
% - gvnames: cell, contains the list of global variables (short names)
% - endog: structure, it contains for each country the corresponding matrix 
%          of endogenous variables
% - endoglist: structure, it contains lists of endogenous variables for each
%              country
% - exog: structure, it contains for each country the corresponding matrix 
%         of weakly exogenous variables
% - exoglist: structure, it contains lists of weakly exogenous variables for 
%             each country
% - dvtype: structure, contains info about type of each variable 
% - wmatrices: structure containing n weight matrices (n: is number of 
%              variable types) 
%--------------------------------------------------------------------------
% OUTPUT:
% - z: struct, contains country series z_it = [domestic_it; foreign_it]
% - x: matrix containing all domestic variables (global vector) 
% - xnames: cell array, contains the corresponding names of variables in x
% - W: struct, contains country link matrices W_i
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************



x=[];  % Initialize global vector xt
xnames=[];


for n=1:cnum
    
    endblock_tmp = endog.(cnames{n});
    endblock = endblock_tmp';
    
    
    exoblock_tmp = exog.(cnames{n});
    exoblock = exoblock_tmp';  
    
    
    % Create country series z (they contain endogenous and foreign
    % variables)
    z.(cnames{n}) = [endblock; exoblock];
    
    x = [ x ; endblock]; %#ok
    xnames=[ xnames; endoglist.(cnames{n})']; %#ok
end

k = rows(x);  % total number of endogenous variables in the system

for i = 1:vnum
    is.(vnames{i}) = zeros(1,k);
end
if not(gvnum==0)
    for i = 1:gvnum
        is.(gvnames{i}) = zeros(1,k);
    end
end


for i=1:k
    for j = 1:vnum
        if strcmp(xnames{i},vnames{j})
            is.(vnames{j})(i) = 1;
        end
    end
    if not(gvnum==0)
        for g = 1:gvnum
            if strcmp(xnames{i},gvnames{g})
                is.(gvnames{g})(i) = 1;
            end
        end
    end
end

% create Link matrices W for each country
sumndom=0;


for n = 1:cnum
    ndom = length(endoglist.(cnames{n})); % # of endogenous variables for each country
    
    % top part of W matrix owing to the x variables 
    wmat_top =[zeros(ndom,sumndom) eye(ndom) zeros(ndom,k-ndom-sumndom)];
    sumndom = sumndom + ndom;
    
    wmat_btm=[];
    wtype = [];
    % bottom part of W matrix owing to the xstar variables
    
        for j=1:fvnum
            fvflag = 0;
            for s=1:length(exoglist.(cnames{n}))
                if strcmp(fvnames(j),exoglist.(cnames{n})(s))
                    fvflag = 1;
                end
            end            
            if fvflag == 1     
                wmat_btm = [wmat_btm; is.(vnames{j})]; %#ok
                wtype = [wtype; dvtype.(fvnames{j})]; %#ok
            end
        end
        pos=0;
        for i=1:cnum 
            for jx = 1:rows(wmat_btm)
               
                wmattype = sprintf('wmat%d',wtype(jx));
                weightsmatrix = wmatrices.(wmattype);
                weight = weightsmatrix(i,n);
               
                wmat_btm(jx,pos+1:pos+length(endoglist.(cnames{i}))) = weight.*wmat_btm(jx,pos+1:pos+length(endoglist.(cnames{i}))); %#ok
            end

            pos = pos + length(endoglist.(cnames{i}));
        end
        
        
        if not(gvnum==0)
            for g=1:gvnum

                gflag = 0;

                for j=1:length(endoglist.(cnames{n}))
                    if strcmp(gvnames(g),endoglist.(cnames{n})(j))
                        gflag = 2;
                    end
                end

                for j=1:length(exoglist.(cnames{n}))
                    if strcmp(gvnames(g),exoglist.(cnames{n})(j))
                        gflag = 1;
                    end
                end
                if gflag == 1
                    wmat_btm = [wmat_btm; is.(gvnames{g})]; %#ok
                end
            end
        end
    
    % let the weights sum to one across the rows
    for i=1:rows(wmat_btm)
        wmat_btm(i,:) = wmat_btm(i,:)./sum(wmat_btm(i,:)); %#ok
    end
    
    % build the link matrices
    W.(cnames{n})=[wmat_top;wmat_btm];
end







