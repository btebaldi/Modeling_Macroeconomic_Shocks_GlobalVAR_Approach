function xvectorized = vec(x)

%**************************************************************************
% PURPOSE: vectorizes a matrix along its columns
%--------------------------------------------------------------------------
% INPUT:
% - x: m x n matrix
%--------------------------------------------------------------------------
% OUTPUT: 
% - xvectorized: (m x n) x 1 array
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
xvectorized = [];

for j=1:cols(x)
    xvectorized = [xvectorized; x(:,j)]; %#ok
end