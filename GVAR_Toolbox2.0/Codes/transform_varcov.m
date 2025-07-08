function varcov = transform_varcov(meth,country_exc,Sigma_u,cnames_s_x)

%**************************************************************************
% PURPOSE: Tranforming the variance-covariance matrix 
%--------------------------------------------------------------------------
% INPUT: 
% - meth: method of computing variance-covariance matrix
%         = 1 do not transform, keep the original estimated varcov matrix
%         = 2 transform the original varcov matrix in a block diagonal one
%         = 3 transform the original varcov matrix in a block diagonal one,
%         with the exception for a particular country specified in
%         country_exc
% - country_exc: exclusion country (short name form) for method 3,it is an
%   empty cell if other methods are chosen
% - Sigma_u: sample varcov matrix
% - cnames_s_x: cell, ordered list of country-labels (short names) of each 
%   endogenous variable in the GVAR (has dimension K times 1) 
%--------------------------------------------------------------------------
% OUTPUT:  
% - varcov: transformed K x K varcov matrix
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014, CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     

if meth == 1 % keep original estimated varcov matrix
    
    varcov = Sigma_u;
    
else  % do transformation
    
    varcov_tmp = Sigma_u;
    for i=1:rows(Sigma_u);
        ci = cnames_s_x{i};
        for j=1:cols(Sigma_u)
            cj = cnames_s_x{j};
                  
            if meth == 2 % block diagonal varcov matrix
                if strcmp(ci,cj)
                else
                    varcov_tmp(i,j) = 0;
                end
            elseif meth == 3 % block diagonal varcov matrix except a specified country
                if strcmp(ci,cj)
                elseif (strcmp(ci,country_exc)) || (strcmp(cj,country_exc))
                else
                    varcov_tmp(i,j) = 0;
                end
            end
            
        end
    end
    varcov = varcov_tmp;
end
