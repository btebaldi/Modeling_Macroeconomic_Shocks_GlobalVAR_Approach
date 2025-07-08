function [ctpcoeffs ctpstd ctpstd_hcw ctpstd_nwc ctptvals ...
         ctptvals_hcw ctptvals_nwc ctplabel]=contmpcoeff(cnum,...
         vnames,endoglist,exoglist,cnames,Psi,std,hcwstd,nwcstd,estcase)

%**************************************************************************
% PURPOSE: retrieves contemporaneous coefficients and corresponding t-values
% from VECMX estimates
%--------------------------------------------------------------------------
% INPUT:
% - cnum: number of countries
% - vnames: cell, list of domestic variables (short names)
% - endoglist: struct, contains lists of endogenous variables for each
% country model
% - exoglist: struct, contains lists of weakly exogenous variables for each
% country model
% - cnames: cell list of countries' names (short names)
% - Psi: struct, contains VECMX* coefficients (excluding the ECM) for each
% country model
% - std: struct, contains standard errors of VECMX* coefficients contained
% in Psi 
% - hcwstd: struct, contains White's Heteroskedasticiy robust standard 
% errors of VECMX* coefficients contained in Psi 
% - nwcstd: struct, contains Newey-West HAC Consistent standard errors
% of VECMX* coefficients contained in Psi 
%--------------------------------------------------------------------------
% OUTPUT: 
% - ctpcoeffs: struct, contains contemporaneous coefficients for each
% country model
% - ctpstd: struct, contains standard errors of ctpcoeffs
% - ctpstd_hcw: struct, contains White's standard errors of ctpcoeffs
% - ctpstd_nwc: struct, contains Newey-West standard errors of ctpcoeffs
% - ctptvals: struct, contains t-values of ctpcoeffs 
% - ctptvals_hcw: struct, contains t-values of ctpcoeffs using White's 
% standard errors
% - ctptvals_nwc: struct, contains t-values of ctpcoeffs using Newey-West 
% standard errors
% - ctplabel: struct, contains variable short names for variable for 
% which the contemporaneous correlation coefficient exists, for each country
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


for n=1:cnum
    
    ndom = length(endoglist.(cnames{n}));
    nfor = length(exoglist.(cnames{n}));
    
    if not(nfor == 0) % skip country models with no foreign-specific variables
        endflag = [];
        exoflag = [];
        for i=1:length(vnames)
            flag=0;
            for j=1:length(endoglist.(cnames{n}))
                if strcmp(vnames(i),endoglist.(cnames{n})(j))
                    flag=1;
                end
            end
            if flag==1
                endflag(i) = 1; %#ok
            else
                endflag(i) = 0; %#ok
            end
        end
        for i=1:length(vnames)
            flag=0;
            for k=1:length(exoglist.(cnames{n}))
                if strcmp(vnames(i),exoglist.(cnames{n})(k))
                    flag=1;
                end
                if flag==1
                    exoflag(i) = 1; %#ok
                else
                    exoflag(i) = 0; %#ok
                end
            end
        end
        
        labflag = zeros(1,length(vnames));
        for i=1:length(vnames)
            if endflag(i)==1 && exoflag(i)== 1
                labflag(i) = 1;
            end
        end
        
        ctplabel.(cnames{n}) = [];
        
        
        for i=1:length(vnames)
            if labflag(i)==1
                ctplabel.(cnames{n}) = [ctplabel.(cnames{n}) vnames(i)];
            end
        end
        
        
        
        inddom = [];
        for h=1:ndom
            for j=1:length(vnames);
                if strcmp(endoglist.(cnames{n})(h),vnames(j))
                    inddom = [inddom h]; %#ok
                end
            end
        end
        
        indfor = [];
        for k=1:nfor
            for j=1:length(vnames);
                if strcmp(exoglist.(cnames{n})(k),vnames(j))
                    indfor = [indfor j]; %#ok
                end
            end
        end
        
        
        
        % extracting coefficients estimates of foreign variables in differences
        pos = cols(Psi.(cnames{n})); % I use a pointer for extracting standard errors
        
        if estcase.(cnames{n}) == 3 || estcase.(cnames{n}) == 4
            inipointer = 2;
            endpointer = 1;
        elseif estcase.(cnames{n}) == 2
            inipointer = 1;
            endpointer = 0;
        end
            
        forcoeffs = Psi.(cnames{n})(:,inipointer:nfor+endpointer);
        stderrors = std.(cnames{n})(:,end-pos+inipointer:end-pos+endpointer+nfor);
        stderrors_hcw = hcwstd.(cnames{n})(:,end-pos+inipointer:end-pos+endpointer+nfor);
        stderrors_nwc = nwcstd.(cnames{n})(:,end-pos+inipointer:end-pos+endpointer+nfor);
        
        ctpcoeffs.(cnames{n}) = [];
        
        ctpstd.(cnames{n}) = [];
        ctpstd_hcw.(cnames{n}) = [];
        ctpstd_nwc.(cnames{n}) = [];
        
        ctptvals.(cnames{n}) = [];
        ctptvals_hcw.(cnames{n}) = [];
        ctptvals_nwc.(cnames{n}) = [];
        
        
        selmat = [];
        for i=1:length(endoglist.(cnames{n}))
            for j=1:length(exoglist.(cnames{n}))
                if strcmp(endoglist.(cnames{n})(i),exoglist.(cnames{n})(j))
                    selmat(i,j) = 1; %#ok
                else
                    selmat(i,j) = 0; %#ok
                end
            end
        end
        
        ctpcoeffs_tmp = selmat.*forcoeffs;
        ctpstd_tmp = selmat.*stderrors;
        ctpstd_hcw_tmp = selmat.*stderrors_hcw;
        ctpstd_nwc_tmp = selmat.*stderrors_nwc;
        
        ctptvals_tmp = ctpcoeffs_tmp./ctpstd_tmp;
        ctptvals_hcw_tmp = ctpcoeffs_tmp./ctpstd_hcw_tmp;
        ctptvals_nwc_tmp = ctpcoeffs_tmp./ctpstd_nwc_tmp;
        
        for i=1:rows(ctpcoeffs_tmp)
            for j=1:cols(ctpcoeffs_tmp)
                if not(ctpcoeffs_tmp(i,j)==0)
                    ctpcoeffs.(cnames{n}) = [ctpcoeffs.(cnames{n}) ctpcoeffs_tmp(i,j)];
                    ctpstd.(cnames{n}) = [ctpstd.(cnames{n}) ctpstd_tmp(i,j)];
                    ctpstd_hcw.(cnames{n}) = [ctpstd_hcw.(cnames{n}) ctpstd_hcw_tmp(i,j)];
                    ctpstd_nwc.(cnames{n}) = [ctpstd_nwc.(cnames{n}) ctpstd_nwc_tmp(i,j)];
                    ctptvals.(cnames{n}) = [ctptvals.(cnames{n}) ctptvals_tmp(i,j)];
                    ctptvals_hcw.(cnames{n}) = [ctptvals_hcw.(cnames{n}) ctptvals_hcw_tmp(i,j)];
                    ctptvals_nwc.(cnames{n}) = [ctptvals_nwc.(cnames{n}) ctptvals_nwc_tmp(i,j)];
                end
            end
        end
    end
end