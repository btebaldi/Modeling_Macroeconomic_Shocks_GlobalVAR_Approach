function [Sigma_zetas Sigma_zeta0 H0_s_t C_s Wy_s beta_s beta_norm_s ncnames ncnames_long ...
         nendoglist nynames ncnames_y ncnames_s_y yorder nyorder] = reorder_GVAR(firstcountries,...
         newordervars,cnamesy,cnamesy_long,Ky,endoglist,...
         H0,C,maxlag,ynames,cnames_s_y,Wy,beta,beta_norm,estcase,Sigma_zeta)
 
%**************************************************************************
% PURPOSE: it changes order of countries and variables in the GVAR
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************     


% Note: the code is general to both SGIRFs and OIRFs
%*****************************************************

nfcnum = length(firstcountries);
cnum = length(cnamesy_long); % note: here the dominant unit model is considered as a country model


n=1;
vlist = {};
for j=1:length(newordervars)
    if not(strcmp(newordervars(j),''))
        vlist = [vlist newordervars(j)]; %#ok
    else
        endoglist_new.(firstcountries{n}) = vlist;
        n=n+1;
        vlist = {};
    end
    
    if j== length(newordervars)
        endoglist_new.(firstcountries{n}) = vlist;
    end
end


% retrieve country long names
%*****************************
firstcountries_long = {};
for j=1:nfcnum
    for n=1:cnum
        if strcmp(firstcountries{j},cnamesy{n})
            firstcountries_long = [firstcountries_long; cnamesy_long(n)]; %#ok
        end
    end
end


% create new cnamesy list
%********************************************
ncnamesy = firstcountries;
ncnamesy_long = firstcountries_long;
for n=1:cnum
    flag = 0;
    for j=1:nfcnum
        if strcmp(cnamesy(n),firstcountries(j))
            flag = 1;
        end
    end
    
    if flag == 0
        ncnamesy = [ncnamesy; cnamesy(n)]; %#ok
        ncnamesy_long = [ncnamesy_long; cnamesy_long(n)]; %#ok
    end
end



% create new endoglist
%********************************************
for n=1:cnum
    flag = 0;
    cpos = 0;
    for j=1:nfcnum
        if strcmp(ncnamesy(n),firstcountries(j))
            flag = 1;
            cpos = j;
        end
    end
    if flag == 1
        nendoglist.(ncnamesy{n}) = endoglist_new.(firstcountries{cpos});
    else
        nendoglist.(ncnamesy{n}) = endoglist.(ncnamesy{n});
    end
end



% create new ynames
%********************************************
nynames=[];
for n=1:cnum
    nynames=[nynames; nendoglist.(ncnamesy{n})']; %#ok
end


% Creating country numbers
%**************************
ncnumbers = zeros(Ky,1);

pos=1;
for n=1:cnum
    
    tmp = length(nendoglist.(ncnamesy{n}));
    ncnumbers(pos:pos+tmp-1) = n;
    
    pos = pos+tmp;
end

% Creating a Ky x 1 cell of country names
%***************************************
ncnames_y = {};
ncnames_s_y = {};

for i=1:Ky
    ncnames_y{i} = ncnamesy_long{ncnumbers(i)}; %#ok
    ncnames_s_y{i} = ncnamesy{ncnumbers(i)}; %#ok
end

ncnames_y = ncnames_y';
ncnames_s_y = ncnames_s_y';


% create a order of variables in the GVAR
%*****************************************
yorder = [1:Ky]'; %#ok
nyorder = zeros(Ky,1);
for i=1:Ky
    vname = nynames(i);
    cname = ncnames_s_y(i);
    for j=1:Ky
        if strcmp(vname,ynames(j)) && strcmp(cname,cnames_s_y(j))
            nyorder(i) = j;
        end
    end
end



% reorder H0 and C matrices
%*************************
H0_t = H0(nyorder,:);
H0_s_t = H0_t(:,nyorder);
for i=1:maxlag
    C_t(:,:,i) = C(nyorder,:,i); %#ok
    C_s(:,:,i) = C_t(:,nyorder,i); %#ok
end

% reorder varcov matrix
%************************
Sigma_zetat = Sigma_zeta(nyorder,:);
Sigma_zetas = Sigma_zetat(:,nyorder);





% reorder beta coefficient vectors and W link matrices
%*****************************************************
if strcmp(firstcountries{1},'du_model')
    k0 = length(nendoglist.du_model);
    vnames_du = ynames(1+end-k0:end)';
    endoglist.du_model = vnames_du;
end

beta_s = beta;
beta_norm_s = beta_norm;


% reorder W link matrices according to the new ordering of variables
for n = 1:cnum 
    Wy_s.(ncnamesy{n}) = Wy.(ncnamesy{n})(:,nyorder);
end

sumk0 = 0;
for nf=1:nfcnum
    k0 = length(nendoglist.(ncnamesy{nf})); % number of endogenous in the first country-models
    sumk0 = sumk0 + k0;
    
    oldlist = endoglist.(ncnamesy{nf});
    newlist = nendoglist.(ncnamesy{nf});
    
    newidx = [];
    for i=1:length(newlist)
        for j=1:length(oldlist)
            if strcmp(oldlist(j),newlist(i))
                newidx = [newidx j]; %#ok
            end
        end
    end
    
    % reorder beta vector
    if isfield(estcase,firstcountries{nf})
        if estcase.(firstcountries{nf}) == 4 || estcase.(firstcountries{nf}) == 2
            betablock = beta_s.(firstcountries{nf})(2:k0+1,:);
            betanormblock = beta_norm_s.(firstcountries{nf})(2:k0+1,:);
        elseif estcase.(firstcountries{nf}) == 3
            betablock = beta_s.(firstcountries{nf})(1:k0,:);
            betanormblock = beta_norm_s.(firstcountries{nf})(1:k0,:);
        end
        
        
        betablock(1:k0,:) = betablock(newidx,:);
        betanormblock(1:k0,:) = betanormblock(newidx,:);
        
        if estcase.(firstcountries{nf}) == 4 || estcase.(firstcountries{nf}) == 2
            beta_s.(firstcountries{nf})(2:k0+1,:) = betablock;
            beta_norm_s.(firstcountries{nf})(2:k0+1,:) = betanormblock;
        elseif estcase.(firstcountries{nf}) == 3
            beta_s.(firstcountries{nf})(1:k0,:) = betablock;
            beta_norm_s.(firstcountries{nf})(1:k0,:) = betanormblock;
        end
    end
    
    % reorder W link matrices 
    Wy_s.(ncnamesy{nf})(1:k0,:) = Wy_s.(ncnamesy{nf})(newidx,:);
end

% recover varcov matrix of selected top countries (which may include the
% dominant unit model)
Sigma_zeta0 = Sigma_zetas(1:sumk0,1:sumk0);

        
% export the new ordering of countries, excluding the dominant unit model if it exists
ncnames = ncnamesy;
ncnames_long = ncnamesy_long;
if strcmp(firstcountries{1},'du_model')
    ncnames = ncnames(2:end);
    ncnames_long = ncnames_long(2:end);
end
    



