function [suplr1 meanlr1 aplr1 suplr2 meanlr2 aplr2 maxobs] = schow(y,x,ccut)


%**************************************************************************
% PURPOSE: do Sequential Chow-tests for breaks. Three tests are constructed
%   (1) Quandt LR = SUP(F)
%   (2) Mean (F)
%   (3) Andrews-Ploberger = ln of mean of (exp(«F))
% and the corresponding heteroskedasticity-robust versions of all three
%--------------------------------------------------------------------------
% INPUT:
% - y = lhv data
% - x = rhv data
% - ccut=endpoints for sequential chow regressions
% - btstrp=indicates whether bootstrap will be conducted
%--------------------------------------------------------------------------
% OUTPUT:
% - suplr1  = sup(f)
% - meanlr1 = mean(f)
% - aplr1   = andrews/ploberger statistic
% - suplr2  = Heteroscedastic Robust Version of sup(f)
% - meanlr2 = Heteroscedastic Robust Version of mean(f)
% - aplr2 = = Heteroscedastic Robust Version of andrews/ploberger statistic
% - maxobs  = observation number corresponding to suplr
%--------------------------------------------------------------------------
% From Gauss code of L. Vanessa Smith. 
% See Dees, di Mauro, Pesaran, Smith (2007).  
%**************************************************************************


nobs=rows(y);
ktrim=floor(ccut*nobs);
n1t=ktrim;
n2t=nobs-ktrim;
lr=zeros(n2t-n1t+1,2);

%         @ full sample  Moment Matrices @
xy=x'*y;
xx=x'*x;
yy=y'*y;


posflag=1;
try
    chol(xx);
catch %#ok
    posflag=0;
end
if posflag==0
xxi = pinv(xx);
else
xxi = eye(rows(xx))/xx;
end


e0=y-x*xxi*xy;
ss0=e0'*e0;
xe=[];
for j=1:cols(x)
    xe_t = x(:,j).*e0;
    xe = [xe xe_t]; %#ok
end 
xxee=xe'*xe;
x1y=x(1:n1t-1,:)'*y(1:n1t-1,:);
x1x1=x(1:n1t-1,:)'*x(1:n1t-1,:);
xxee1=xe(1:n1t-1,:)'*xe(1:n1t-1,:);
yy1=y(1:n1t-1,:)'*y(1:n1t-1,:);


i=n1t;
while i<= n2t;

    x1y=x1y+x(i,:)'*y(i,:);
    x1x1=x1x1+x(i,:)'*x(i,:);
    xxee1=xxee1+xe(i,:)'*xe(i,:);
    yy1=yy1+y(i,:)'*y(i,:);
    x2x2=xx-x1x1;
    x2y=xy-x1y;
    yy2=yy-yy1;
    xxee2=xxee-xxee1;

    
    posflag=1;
    try
        chol(x1x1);
    catch %#ok
        posflag=0;
    end
    if posflag==0
        x1x1i = pinv(x1x1);
    else
        x1x1i = eye(rows(x1x1))/x1x1;
    end
    
    posflag=1;
    try
        chol(x2x2);
    catch %#ok
        posflag=0;
    end
    if posflag==0
        x2x2i = pinv(x2x2);
    else
        x2x2i = eye(rows(x2x2))/x2x2;
    end

    b1=x1x1i*x1y;
    b2=x2x2i*x2y;

    ssb=(yy1-x1y'*b1) + (yy2-x2y'*b2);

    lr(i-n1t+1,1)=(nobs-2*cols(x))*((ss0-ssb)/ssb);

    %           @ -- Hetero- Robust F-test -- @

    v=x1x1i*(xxee1)*x1x1i + x2x2i*(xxee2)*x2x2i;

    
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
    
    lr(i-n1t+1,2)=(b1-b2)'*v1*(b1-b2);

    i=i+1;
end

lr1=lr(:,1);
lr2=lr(:,2);

% @ -- Non Robust -- @
suplr1=maxc(lr1);
maxobs=n1t+maxindc(lr1)-1;
meanlr1=mean(lr1);
aplr1=log(mean(exp(0.5*lr1)));

%   @ -- Hetero Robust -- @
suplr2=maxc(lr2);
meanlr2=mean(lr2);
aplr2=log(mean(exp(0.5*lr2)));
