function xforc = forecast_GVAR(fhorz,sposfirst,sposlast,K,x,maxlag,C,delta_0,delta_1,...
                            trIR,virlb,xnames,lb_flag,lb_restr)

%**************************************************************************
% PURPOSE: It computes ex-ante forecasts of the GVAR model subject to the 
% lower bound restriction of 0.25 on the per annum interest rates (short  
% and long) of all countries, should such interest rates be present in the 
% GVAR model and identified by the user via the interface file. In addition
% it allows the user to input additional restrictions on additional country-
% specific variables (these do not have to be imposed on all countries). 
%**************************************************************************     

x_lb=-88888888*ones(size(xnames,1),1);

if isempty(virlb)==0;
    
    indx=zeros(size(xnames,1),1);
    
    for i=1:size(virlb,1);
        indx=indx+strcmp(xnames,virlb(i));
    end
    
    idvec=find(indx==1);
    
    x_lb(idvec)=trIR*ones(sum(indx),1); %#ok
    
end


if lb_flag==1;
    
    vaux=find(~isnan(lb_restr));
    
    if isempty(vaux)==0;
        x_lb(vaux)=lb_restr(vaux);
    end
end


indxyy = sposlast-sposfirst+1;
sumki = K;
xt_tr = x;
pmx = maxlag;

GBefwm = [];
for i=1:pmx
    GBefwm = [GBefwm C(:,:,i)]; %#ok
end

xt=xt_tr(:,(indxyy-(pmx-1)):cols(xt_tr)); % actual out of sample period,
% including maxlag additional actual values to initialise the forecasting

xlags = xt(:,1:pmx);

rsq_tmp = cols(xlags)*ones(cols(xlags),1);
rsq = rsq_tmp - (0:1:cols(xlags)-1)';

xlagsr = xlags(:,rsq); % columns reversed, for convenience later

xforc = -9999*ones(sumki,fhorz);

conest = delta_0; % estimated GVAR intercept
trndest = delta_1; % estimated GVAR trend

jj=1;
while jj<=fhorz
    
    sme = zeros(sumki,1);
    jwm = 1;
    while jwm <= pmx
        sme = sme + GBefwm(:,(jwm-1)*sumki+1:(jwm-1)*sumki+sumki)*xlagsr(:,jwm);
        jwm=jwm+1;
    end
    
    trnd = 0:1:cols(xt_tr)-1;
    
    xforcjj = conest+trndest*(trnd(indxyy)+jj)+sme;
    
    compind=find(xforcjj < x_lb);
    if isempty(compind)~=1;
        xforcjj(compind)=x_lb(compind);
    end
    
    xforc(:,jj) = xforcjj;
    
    if pmx>1
        xlagsr=[xforc(:,jj) xlagsr(:,1:pmx-1)];
    else
        xlagsr=xforc(:,jj);
    end
    
    jj=jj+1;
end