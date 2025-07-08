function PPres = pprofile(PHI,Sigma_u,G,W,beta,N,cnames,estcase)

%**************************************************************************
% PURPOSE: Computing persistence profiles
%--------------------------------------------------------------------------
% INPUT:
% - PHI: three-dim matrix (K x K x (N+1)) containing dynamic multipliers of
% the GVAR
% - Sigma_u: K x K varcov matrix of the GVAR
% - G: K x K matrix containing impact coefficients of the GVAR
% - W: struct, contains country link matrices W_i 
% - beta: struct, contains cointegrating vector coefficients estimates
% - N: forecast horizon
% - cnames: list of countries' names (short names)
% - cnum: number of countries
% - estcase: struct, containstreatment of deterministic components in the 
% estimation: if =4,it is case IV (restricted trend in cointegration space,
% unrestricted intercept in levels); if =3, it is case III (no trend in
% cointegration space, unrestricted intercept in levels); if = 2, it is case
% II (no trend, restricted intercept in cointegration space). Cases are from
% MacKinnon, Haug and Michelis (1999)
%--------------------------------------------------------------------------
% OUTPUT:
% - PPres: struct, contains persistence profiles for each cointegrating
% vector of each country in the GVAR
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************

cnum = length(cnames);



PPres=[];

invG = eye(rows(G))/G;
invGSigma = G\Sigma_u;

for n=1:cnum
    if not(isempty(beta.(cnames{n})))
        PPres_idv = [];
        
        Wm = W.(cnames{n});
        
        if estcase.(cnames{n}) == 4 || estcase.(cnames{n}) == 2
            betablock = beta.(cnames{n})(2:end,:);
        elseif estcase.(cnames{n}) == 3
            betablock = beta.(cnames{n});
        end
        
        den = betablock'*Wm*invGSigma*invG'*Wm'*betablock;
        
        for i=1:N+1
            num = betablock'*Wm*PHI(:,:,i)*(invGSigma)*invG'*PHI(:,:,i)'*Wm'*betablock;
            
            PPres_temp = num./den;
            
            PPres_idv = [PPres_idv diag(PPres_temp)]; %#ok
        end
        
        PPres = [PPres; PPres_idv]; %#ok
    end
end
