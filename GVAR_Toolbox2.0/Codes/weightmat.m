function weights = weightmat(trdmtrx)

%**************************************************************************
% PURPOSE: create a matrix of weights 
%--------------------------------------------------------------------------
% INPUT: 
% - trdmtrx: a k x k matrix of levels (e.g. imports/exports)
%--------------------------------------------------------------------------
% OUTPUT: 
% - trdwghts: a k x k matrix of trade weights
%--------------------------------------------------------------------------
% NOTE: the matrix is created such that columns, but not rows, sum to one
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
weights=zeros(rows(trdmtrx),cols(trdmtrx));

for i=1:rows(trdmtrx)
    for j=1:cols(trdmtrx)
        if i~=j
            weights(i,j)= trdmtrx(i,j)/sum(trdmtrx(:,j));
        else 
        end
    end
end
