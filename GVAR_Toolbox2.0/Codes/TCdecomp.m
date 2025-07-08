function [xp xc x_tilde xp_dt xp_st a] = TCdecomp(TC_RestrictionFlag,F,...
                                         eta,maxlag,x,x_resIndicator,cond_TC_noDU,cond_TC_DU,dumodel_flag)
                           
%**************************************************************************
% PURPOSE: It decomposes the variables in the GVAR into permanent and
%       transitory components. 
%--------------------------------------------------------------------------
% INPUT:
% - TC_RestrictionFlag:
% - F: three-dim matrix (K x K x maxlag) containing estimated 
% coefficients of the GVAR (i.e. of reduced form representation)
% - eta: K x nobs matrix containing residuals of the reduced-form GVAR 
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
% - x: matrix containing all domestic variables (global vector) 
% - x_resIndicator: vector of trend restrictions
%--------------------------------------------------------------------------
% OUTPUT:
% - X_permanent: permanent component  
% - X_transitory: transitory/cyclical component 
% - X_tilda: X-X_permanent (should be equal to the cyclical component)
%--------------------------------------------------------------------------
% Meiling He and L. Vanessa Smith, September 2011
% CFAP, Judge Business School, University of Cambridge
% mh590@cam.ac.uk and lvs21@cam.ac.uk
%**************************************************************************

% disp('---- * Performing the Trend/Cycle Decomposition of the GVAR * ----');
% %1) Calculate C(1)
% disp('- Calculating C(1)');

[k T]=size(x);
%k = rows(x);
%T = cols(x);
sumLimit = 1000;
C = zeros([k,k,sumLimit]);
C(:,:,1) = eye(k);                   %C(0)
C(:,:,2) = - (C(:,:,1) - F(:,:,1));         % C(1)

for j = 3:sumLimit  % C(2) to c(1000)
    for n = 1:maxlag
        if j-n>0
            c_t = C(:,:,(j-n))*F(:,:,n);
            C(:,:,j)=C(:,:,j) + c_t; 
        end
    end
end

C1 = sum(C,3); % C(1)

%2) Calculate the Permanent-stochastic component, Xp(st)
% disp('- Calculating the permanent-stochastic component, Xp(st)')
xp_st = zeros([k,T-maxlag]);
for t = 1:(T-maxlag)
    sum_error = sum(eta(:,1:t),2);
    xp_st(:,t) = C1*sum_error;
end

%3) Calculate vt
% disp('- Calculating vt')
xt = x(:,maxlag+1:T);
vt = xt - xp_st;

%4)OLS estimation 
%4.1)Creat trend
trend_t = [0:1:T]';
trend = (trimr(trend_t,maxlag,1));

%4.2)creat indicator for restricted variables
    
%4.3)OLS without/with restriction 

if TC_RestrictionFlag==1 %OLS with restriction
    
    x_lineIndex = [1:k]';
    x_noRes_lineIndex = x_lineIndex(~x_resIndicator);
    x_Res_lineIndex =  x_lineIndex( x_resIndicator);
   
    x_OLS1 = [ones(T-maxlag,1) trend];
    a1 = (x_OLS1)\(vt(x_noRes_lineIndex,:))';
    
    x_OLS2 = [ones(T-maxlag,1)];
    a2 = (x_OLS2)\(vt(x_Res_lineIndex,:))';
   
    a(:,x_noRes_lineIndex) = a1;
    a(1,x_Res_lineIndex) = a2;
    a(2,x_Res_lineIndex) = 0;

else %without restriction
%     disp('- OLS without restriction ')
       
  if dumodel_flag==0           
     if cond_TC_noDU==1
     x_OLS = [ones(T-maxlag,1)];  
     else    
     x_OLS = [ones(T-maxlag,1) trend];
     end       
  elseif dumodel_flag==1          
     if cond_TC_DU==1 
     x_OLS = [ones(T-maxlag,1)];   
     else    
     x_OLS = [ones(T-maxlag,1) trend];
     end           
  end
  
   a = x_OLS\vt';

end 

%4.4) Calculating xc, xp, and x_tilda
% disp('- Calculating xc, xp, and x_tilda ')

if dumodel_flag==0
  if cond_TC_noDU==1
  xp_dt = a(1,:)'*(ones(T-maxlag,1))';    
  else    
  xp_dt = a(1,:)'*(ones(T-maxlag,1))' + a(2,:)'*trend';
  end
elseif dumodel_flag==1         
  if cond_TC_DU==1 
  xp_dt = a(1,:)'*(ones(T-maxlag,1))';    
  else    
  xp_dt = a(1,:)'*(ones(T-maxlag,1))' + a(2,:)'*trend';
  end    
end    

xp = xp_st + xp_dt;

xc = vt - xp_dt;

x_tilde = xt - xp; 

% disp('- Comparing xc and x_tilda ')
comp_t = xc-x_tilde;

comp_trend = sum(comp_t);

comp_sum = sum(comp_trend );

if (comp_sum>0.0000000001)
    disp('');
    disp('- The cyclical component and the deviations from the permanent component are not identical.');
    disp('  Please contact the GVAR helpdesk: gvar.helpdesk@gmail.com ');
    error('- End of program. ')
    
else
%     disp('')
%     disp('-- No significant difference between x_tida and xc ----')
    
end

% disp('---- * TC Decomposition completed * ----')











