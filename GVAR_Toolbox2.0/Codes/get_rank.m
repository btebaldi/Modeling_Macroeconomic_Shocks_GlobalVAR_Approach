
function [rank trace_critvals maxeig_critvals] = get_rank(cnum,cnames,...
    endog,exog,trace,trace_crit95_c2,trace_crit95_c3,trace_crit95_c4,maxeig,...
    maxeig_crit95_c2,maxeig_crit95_c3,maxeig_crit95_c4,estcase)

%**************************************************************************
% PURPOSE: It yields the number of cointegrating relations for each country 
%--------------------------------------------------------------------------
% INPUT:
% - cnum: number of country models
% - cnames: cell, list of countries' names (short names)
% - endog: struct, it contains for each country the corresponding matrix 
% of endogenous variables
% - exog: struct, it contains for each country the corresponding matrix 
% of weakly exogenous variables
% - trace: struct, it contains the trace statistics for each country model
% - trace_crit95_c2: matrix containing critical values (95% confidence
% level) for trace test, for case II (no trend, unrestricted intercept in cointegration space)
% - trace_crit95_c3: matrix containing critical values (95% confidence
% level) for trace test, for case III (no trend in coint space,
% unrestricted intercepts in levels
% - trace_crit95_c4: matrix containing critical values (95% confidence
% level) for trace test, for case IV (restricted trend in coint space,
% unrestricted intercepts in levels)
% - maxeig: struct, it contains the maximum eigenvalue statistics for each
% country model
% - maxeig_crit95_c2: matrix containing critical values (95% confidence
% level) for maximum eigenvalue test, for case II (no trend, unrestricted intercept in
% cointegration space)
% - maxeig_crit95_c3: matrix containing critical values (95% confidence
% level) for maximum eigenvalue test, for case III (no trend in coint space,
% unrestricted intercepts in levels)
% - maxeig_crit95_c4: matrix containing critical values (95% confidence
% level) for maximum eigenvalue test, for case IV (restricted trend in 
% coint space, unrestricted intercepts in levels)
% - estcase: struct, it contains for each country a value related to the
% treatment of deterministic components in the estimation: 
% if =4, it is case IV (restricted trend in cointegration space,
% unrestricted intercept in levels); if =3, it is case III (no trend in
% cointegration space, unrestricted intercept in levels). Cases are from
% MacKinnon, Haug and Michelis (1999)
%--------------------------------------------------------------------------
% OUTPUT:
% - rank: struct, it contains ranks of country models
% - trace_critvals: struct, it contains trace statistic's critical values
% for each country model
% - maxeig_critvals: struct, it contains maximum eigenvalue statistic's
% critical values for each country model
%--------------------------------------------------------------------------
% NOTES: 
% - by default, ranks are determined using Trace statistic 
% - critical values are at 95% significance level
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


% Settings:
% choose stat: 0 for trace, 1 for max eigenvalue
stat = 0; 
% (fixed) level of significance: 95%


for n = 1: cnum
    rank.(cnames{n})=0;
    endnum = cols(endog.(cnames{n}));
    exonum = cols(exog.(cnames{n}));
    
    if estcase.(cnames{n}) == 4
        
        critvals_tmp = trace_crit95_c4(1:endnum, exonum+1);
        trace_critvals.(cnames{n}) = flipud(critvals_tmp);
        
        critvals_tmp = maxeig_crit95_c4(1:endnum, exonum+1);
        maxeig_critvals.(cnames{n}) = flipud(critvals_tmp);
        
    elseif estcase.(cnames{n}) == 3
        
        critvals_tmp = trace_crit95_c3(1:endnum, exonum+1);
        trace_critvals.(cnames{n}) = flipud(critvals_tmp);
        
        critvals_tmp = maxeig_crit95_c3(1:endnum, exonum+1);
        maxeig_critvals.(cnames{n}) = flipud(critvals_tmp);
        
    elseif estcase.(cnames{n}) == 2
        
        critvals_tmp = trace_crit95_c2(1:endnum, exonum+1);
        trace_critvals.(cnames{n}) = flipud(critvals_tmp);
        
        critvals_tmp = maxeig_crit95_c2(1:endnum, exonum+1);
        maxeig_critvals.(cnames{n}) = flipud(critvals_tmp);
        
    end
    
    if stat == 0  % Trace statistic chosen
        for i=1:endnum
            if trace.(cnames{n})(i) > trace_critvals.(cnames{n})(i)
                rank.(cnames{n})=rank.(cnames{n})+1;
            else
                break
            end
        end
    elseif stat == 1  % Maximum eigenvalue statistic chosen
        for i=1:endnum
            if maxeig.(cnames{n})(i) > maxeig_critvals.(cnames{n})(i)
                rank.(cnames{n})=rank.(cnames{n})+1;
            else
                break
            end
        end
    end
end
