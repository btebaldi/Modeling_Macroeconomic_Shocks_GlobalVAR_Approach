function res = maxc(x)

%**************************************************************************
% PURPOSE: Given a matrix, it yields maximum values over columns, and
% resulting output is a column vector (rather than a row vector as the
% MatLab routine max()
%--------------------------------------------------------------------------
% INPUT:
% - x: m x n matrix
%--------------------------------------------------------------------------
% OUTPUT:
% - res: n x 1 matrix of maximum values over each column of x
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

res_t = max(x);

res = res_t';