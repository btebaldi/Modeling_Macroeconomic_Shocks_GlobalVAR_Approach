
function [seNW q] = neweywest(X,res)

%**************************************************************************
% PURPOSE: computes Newey-West Heteroskedasticity and Autocorrelation
% Consistent standard errors
%--------------------------------------------------------------------------
% From Gauss code of L. Vanessa Smith. 
% See Dees, di Mauro, Pesaran, Smith (2007).  
%**************************************************************************



s = cols(X);
T = rows(res);
q = floor(4*(T/100)^(2/9));  % following Newey-West suggestion

for v=1:q
          cw=zeros(s,s);
          iw=v+1;
          while iw<=rows(X);
          qi=X(iw,:)'*res(iw,1)*res(iw-v,1)*X(iw-v,:);
          cw=cw+qi;
          iw=iw+1;
          end
   
cw_new=(1-(v/(q+1)))*(cw+cw');

if v==1;
cw_new1=cw_new;
elseif v>1;
cw_new1=cw_new1+cw_new;
end

end

lalagsbda= diag(diag(res*res'));

invXX = eye(rows(X'*X))/(X'*X);

varcovNW= (T/(T-s))*(invXX*X'*lalagsbda*X*invXX+invXX*cw_new1*invXX); 

seNW=sqrt(diag(varcovNW));