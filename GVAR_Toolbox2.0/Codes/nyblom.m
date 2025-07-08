
function [lm rlm] =nyblom(y,x2)

%**************************************************************************
% PURPOSE: computes the Nyblom (1989) test for parameter constancy against
% nonstationary alternatives (here lm) and also its
% heteroskedasticity-robust version (here rlm)
%--------------------------------------------------------------------------
% From Gauss code of L. Vanessa Smith. 
% See Dees, di Mauro, Pesaran, Smith (2007).  
%**************************************************************************


yy = y;
zz = x2;
k = cols(x2);


% check if zz'*zz is posdef: if so avoid doing pseudoinverse
posflag=1;
try
    chol(zz'*zz);
catch %#ok
    posflag=0;
end

if posflag==0
mzinv = pinv(zz'*zz);
else
mzinv = eye(rows(zz'*zz))/(zz'*zz);
end


e=yy-zz*(mzinv*(zz'*yy));
seesq=(e'*e)/(rows(e)-k);
ex=[];
for j=1:cols(zz)
    ex_t = zz(:,j).*e;
    ex = [ex ex_t]; %#ok
end
exs=cumsum(ex);
 
v=seesq*(zz'*zz);
rv=ex'*ex;


posflag=1;
try
    chol(v);
catch %#ok
    posflag=0;
end
if posflag==0
v1 = pinv(v);
else
v1 = eye(rows(v))/v;
end

posflag=1;
try
    chol(rv);
catch %#ok
    posflag=0;
end
if posflag==0
v2 = pinv(rv);
else
v2 = eye(rows(rv))/rv;
end


lm=(v1*(exs'*exs))/rows(e);
dglm=diag(lm);
lm=sum(dglm);

rlm=(v2*(exs'*exs))/rows(e);
dgrlm=diag(rlm);
rlm=sum(dgrlm);
