
function [m1 m2] = kraplob(y,x)

%**************************************************************************
% PURPOSE: computes the Ploberger and Kramer's (1992) maximal OLS
% cumulative sum (CUSUM) statistic, here callede m1, and its mean-square 
% variant, here called m2.
%--------------------------------------------------------------------------
% From Gauss code of L. Vanessa Smith. 
% See Dees, di Mauro, Pesaran, Smith (2007).  
%**************************************************************************

% Galesi:
% check if x'*x is posdef: if so avoid doing pseudoinverse
posflag=1;
try
    chol(x'*x);
catch %#ok
    posflag=0;
end

if posflag==0
invxx = pinv(x'*x);
else
invxx = eye(rows(x'*x))/(x'*x);
end

ehat=y - x*invxx*(x'*y);
vcv=ehat'*ehat/(rows(y)-cols(x));
t1=cumsum(ehat/sqrt(rows(ehat)*vcv)); 
m1=max(abs(t1));
t2=t1.*t1; 
m2=mean(t2);

