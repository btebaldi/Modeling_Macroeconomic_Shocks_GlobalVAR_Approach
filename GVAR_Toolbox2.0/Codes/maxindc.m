function res = maxindc(x)

%**************************************************************************
% PURPOSE: Given a matrix, it yields positions of maximum values over 
% columns 
%--------------------------------------------------------------------------
% INPUT:
% - x: m x n matrix
%--------------------------------------------------------------------------
% OUTPUT:
% - res: n x 1 matrix of positions of maximum values over each column of x
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

[junk indres_t] = max(x); %#ok

res = indres_t';