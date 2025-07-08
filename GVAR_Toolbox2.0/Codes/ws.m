function wsstat = ws(x,d,p)

%**************************************************************************
% PURPOSE: Performs Weighted Symmetric Dickey-Fuller Test for Unit Root
% Presence
%--------------------------------------------------------------------------
% INPUT:
% - y: a nobs x 1 time series
% - d: indicates the type of detrending to perform before undertaking the
%      WS unit root tests
%      if d=0, regression of x under an intercept 
%      if d=1, regression of x under an intercept and a linear trend
% - p: # of lagged changes to include in the fitted regression 
%--------------------------------------------------------------------------
% OUTPUT: 
% - wsstat: WS statistic
%--------------------------------------------------------------------------
% From Gauss code of L. Vanessa Smith. 
% See Dees, di Mauro, Pesaran, Smith (2007).  
%**************************************************************************


x = detrend(x,d);

if d==0
    c=2;
elseif d==1
    c=3;
end

T = rows(x);

sum1a=zeros(p+1,p+1);
sum1b=zeros(p+1,1);

Dx = trimr(x-lagm(x),1,0);
k=0;
zbs = trimr(lagm(x),1,0);

if p>0
    while k < p    
        k = k+1;
        zbs = [zbs lagm(Dx,k)]; %#ok
    end
end

zbs = trimr(zbs,k,0);


w1=-9999*ones(rows(x),1);

s1=1;
while s1 <= T
    if (1<=s1 && s1<=p+1)
        w1(s1)=0;
    elseif ((p+1)<s1) && (s1<=T-p)
        w1(s1) = (s1-p-1)/(T-2*p);
    elseif (T-p<s1<=T)
        w1(s1) = 1;
    end
    s1 = s1+1;
end

    w1tr  = trimr(w1,k+1,0);
    xbs =trimr(x,k+1,0);
    
i1=1;
while i1 <= rows(zbs)
    sum1a=sum1a+w1tr(i1)*(zbs(i1,:)'*zbs(i1,:));
    sum1b=sum1b+w1tr(i1)*zbs(i1,:)'*xbs(i1);
   
    i1 = i1+1;
end

w2=trimr(w1,1,p);

sum2a=zeros(p+1,p+1);
sum2b=zeros(p+1,1);   
       
     k       = 0;
  zfs=trimr(x,1,0);
  
if p>0
    zfs = trimr(zfs,0,p);
    while k<p              %%%%%%
        k=k+1;
        xksin1=trimr(x,k+1,p-k);
        xk=trimr(x,k,p-(k-1));
        
        zfs = [zfs -(xksin1-xk)]; %#ok
    end
end

depf     = trimr(x,0,k+1);

i2=1;
while i2 <= rows(zfs)

sum2a=sum2a+(1-w2(i2))*(zfs(i2,:)'*zfs(i2,:));
sum2b=sum2b+(1-w2(i2))*zfs(i2,:)'*x(i2);

i2=i2+1;
end


As=sum1a+sum2a;

bs=sum1b+sum2b;

theta=As\bs;
        

sum1q=0;
q1=1;
while q1<=rows(zbs)

if p>0

indb=[xbs zbs];

sum1q=sum1q+w1tr(q1)*([1 (-1)*theta']*indb(q1,:)')^2;

elseif p==0

indb=[xbs zbs(:,1)];

sum1q=sum1q+w1tr(q1)*([1 (-1)*theta(1)]*indb(q1,:)')^2;

end

q1=q1+1;
end       



sum2q=0;
q2=1;
while q2<=rows(zfs)

if p>0;

indf=[depf zfs];

sum2q=sum2q+(1-w2(q2))*([1 (-1)*theta(1) ((-1)*theta(2:rows(theta)))']*indf(q2,:)')^2;

elseif p==0;

indf=[depf zfs(:,1)];

sum2q=sum2q+(1-w2(q2))*([1 (-1)*theta(1)]*indf(q2,:)')^2;

end

q2=q2+1;
end


Q=sum1q+sum2q;


Ainv=eye(rows(As))/As;

estvar=(Q/(T-p-c))*Ainv(1,1);



wsstat=(theta(1,1)-1)/sqrt(estvar);

