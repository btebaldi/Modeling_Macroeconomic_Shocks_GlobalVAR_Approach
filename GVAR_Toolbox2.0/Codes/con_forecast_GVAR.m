function [mu_s] = con_forecast_GVAR(maxlag,sposfirst,sposlast,...
    delta_0,delta_1,x,K,Sigma_eta,C,...
    con_fhorz,con_fhorz_restr,...
    con_forc_restr_mtx,con_forc_restr_x,xnames)



%**************************************************************************
% PURPOSE: calculates the conditional forecasts for the GVAR
%--------------------------------------------------------------------------
% L. Vanessa Smith, Nov 2012
% CFAP, University of Cambridge
% lvs21@cam.ac.uk
%**************************************************************************
%

d_Hbar_mtx = con_forc_restr_mtx;%

nbrestr = sum(con_forc_restr_x);

Psi = zeros(nbrestr,K);
Psi(:,con_forc_restr_x) = eye(nbrestr);

mu_s = -9999*ones(K,con_fhorz);

covmtx = zeros(maxlag*K);

covmtx(1:K,1:K) = Sigma_eta;

Psi_kr = kron(eye(con_fhorz_restr),Psi);

g_Hbar_tilda_mtx = -9999*ones(nbrestr,con_fhorz_restr);

omega_Hbar_tilda = -9999*ones(K*con_fhorz_restr);

Cf  = zeros(maxlag*K); %Cf is the companion form of C
for i=1:maxlag
    Cf(1:K,(1+K*(i-1)):(K+K*(i-1))) = C(:,:,i);
end
for i = 1:maxlag-1
    Cf((1+K*i):(K+K*i),(1+K*(i-1)):(K+K*(i-1))) = eye(K);
end 
       
sm = zeros(rows(Cf),cols(Cf));

Cfi = eye(rows(Cf));


trIR=[];
virlb=[];
lb_flag=0;
lb_restr=[];

mu = forecast_GVAR(con_fhorz_restr,sposfirst,sposlast,K,x,maxlag,C,delta_0,delta_1,...
         trIR,virlb,xnames,lb_flag,lb_restr);

     
for i = 1:con_fhorz_restr
    ss1=(i-1)*K+1;
    ss2=(i-1)*K+K;
    
    sm=sm+Cfi*covmtx*Cfi';
    Cfi = Cfi*Cf;
    g_Hbar_tilda_mtx(:,i) = d_Hbar_mtx(:,i)-Psi*mu(:,i);
    
    Cfnew = Cf;
    
    for j = 1:con_fhorz_restr
        if i==j
            omega_Hbar_tilda(ss1:ss2,ss1:ss2)=sm(1:(rows(Cf)/maxlag),1:(cols(Cf)/maxlag));
        elseif i<j
            ss3=(j-1)*K+1;
            ss4=(j-1)*K+K;
            sm1new = sm*Cfnew';
            sm2new = Cfnew*sm;
            
            omega_Hbar_tilda(ss1:ss2,ss3:ss4)=sm1new(1:(rows(Cf)/maxlag),1:(cols(Cf)/maxlag));
            omega_Hbar_tilda(ss3:ss4,ss1:ss2)=sm2new(1:(rows(Cf)/maxlag),1:(cols(Cf)/maxlag));
            
            Cfnew = Cfnew*Cf;
        end
    end
end

for h=1:con_fhorz
    sg1 = (h-1)*K + 1;
    sg2 = (h-1)*K + K;
    
    omega_Hbar_hL = omega_Hbar_tilda(sg1:sg2,:);
    
    mu_s(:,h) = mu(:,h) + omega_Hbar_hL*Psi_kr'*inv(Psi_kr*omega_Hbar_tilda*Psi_kr')*vec(g_Hbar_tilda_mtx); %#ok
    
end 



