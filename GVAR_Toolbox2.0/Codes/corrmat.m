function cmat = corrmat(x)


%**************************************************************************
% PURPOSE: computes pairwise linear correlations
%--------------------------------------------------------------------------
% INPUT: 
% - x: T x P matrix (say T is the number of temporal observations, P the
% number of variables)
%--------------------------------------------------------------------------
% OUTPUT:
% - cmat: P x P matrix containing the pairwise linear
%     correlation coefficient between each pair of columns in the T x P
%     matrix x.
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************


T = rows(x);
k = cols(x);

cmat = zeros(k);

for i=1:k
    for j=1:k
        p1 = x(:,i) - ones(T,1)*mean(x(:,i));
        p2 = x(:,j) - ones(T,1)*mean(x(:,j));

        var1 = p1'*p1;
        var2 = p2'*p2;

        cov12 = p1'*p2;

        corr12 = cov12/(sqrt(var1)*sqrt(var2));
        
        cmat(i,j) = corr12;
        cmat(j,i) = corr12;
    end
end