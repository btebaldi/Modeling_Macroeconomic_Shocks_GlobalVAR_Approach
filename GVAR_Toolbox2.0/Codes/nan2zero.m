function y = nan2zero(x)

%**************************************************************************
% PURPOSE: Substituting missing values with zeros 
%--------------------------------------------------------------------------
% INPUT: x: a m x n matrix with some NaNs
%--------------------------------------------------------------------------
% OUTPUT: y: a m x n matrix with zeros in place of NaNs
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

y=x;

for i=1:rows(x)
    for j=1:cols(x)
        if isnan(x(i,j)) == 1
            y(i,j)=0;
        end
    end
end

