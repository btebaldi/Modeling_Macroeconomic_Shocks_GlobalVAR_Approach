% for i=1:137
%     modelName = strcat('mdl', num2str(i));
% disp(modelName);
% 
%     X=a(:,i);
%     Y=b(:,i);
%     Models.(modelName)=fitlm(X,Y);
%     
% end

clear all;

load('Teste.mat')

clc

% Y=[adm desl adms desls dlipca lpim_BR lsel];
Y=[adm desl adms desls];

[Johansen_h, Johansen_pValue, Johansen_stat, Johansen_cValue, Johansen_mles] = jcitest(Y, 'lags', 1);

mdl = fitlm([desls], adms, ...
    'ResponseVar','Admit*',...
    'PredictorVars',{'Deslig*'},'Intercept',false)

beta = table2array(mdl.Coefficients(1,1));

new = adms-1.0427*desls;

fitlm([desl new], adm, ...
    'ResponseVar','Admitidos',...
    'PredictorVars',{'Deslidagos','coint'},'Intercept',false)

