function [newtradematrix newcnames_long newcnames newcnum] = update_matrix(cnum,cnames,cnames_long,name,name_long,list,tradematrix)


%**************************************************************************
% PURPOSE: resize the weights matrix according to the created region 
%--------------------------------------------------------------------------
% INPUT:
% - cnum: number of countries
% - cnames: list of countries' names (short names)
% - cnames_long: list of countries' names (long names)
% - name: name of the region that has been constructed (short name)
% - name_long: name of the region that has been constructed (long name)
% - list: list of countries belonging to region 'name'
% - tradematrix: cnum x cnum matrix of levels (e.g. trade flows)
%--------------------------------------------------------------------------
% OUTPUT:
% - newtradematrix: newcnum x newcnum matrix of levels (e.g. trade flows)
% - newcnames_long: updated list of countries' names (long names)
% - newcnames: updated list of countries' names (short names)
% - newcnum: updated number of countries
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
regnum = length(list);

regind = zeros(regnum,1);
allind = zeros(cnum,1);
restind = zeros(cnum-regnum,1);

k=1; h=1; 
for n=1:cnum
    if strcmp(cnames{n},list{k})
        regind(k) = n;
        allind(n) = 1;
        if k==regnum
        else
        k=k+1;
        end
    else
        restind(h) = n;
        if h==(cnum-regnum)
        else
            h=h+1;
        end
    end
end



newmat=zeros(cnum,cnum);
for i=1:rows(tradematrix)
    flagi=0;
    for k=1:length(regind)
        if i==regind(1)
            flagi=1;
        elseif ((i==regind(k) && not(i==regind(1))))
            flagi=2;
        end
    end

    if flagi ==0
        newmat(i,1)=tradematrix(i,1);
        for j=1:cols(tradematrix)
            flagj=0;
            for k=1:length(regind)
                if j==regind(1)
                    flagj=1;
                elseif ((j==regind(k) && not(j==regind(1))))
                    flagj=2;
                end
            end

            if flagj == 0
                newmat(i,j) = tradematrix(i,j);
            elseif flagj ==1
                acc = 0;
                for k=1:length(regind)
                    acc=acc+tradematrix(i,regind(k));
                end
                newmat(i,j)= acc;
            elseif flagj == 2
                newmat(i,j) = NaN;
            end
        end
    elseif flagi==1
        acc=0;
        for k=1:length(regind)
            acc=acc+tradematrix(regind(k),1);
        end
        newmat(i,1) = acc;

        for j=1:cols(tradematrix)
            flagj=0;
            for k=1:length(regind)
                if j==regind(1)
                    flagj=1;
                elseif ((j==regind(k) && not(j==regind(1))))
                    flagj=2;
                end
            end

            if flagj==0
                acc=0;
                for k=1:length(regind)
                    acc=acc+tradematrix(regind(k),j);
                end
                newmat(i,j) = acc;

            elseif flagj==1
                newmat(i,j) = tradematrix(i,j);
            elseif flagj==2
                newmat(i,j) = NaN;
            end
        end
    elseif flagi==2
        newmat(i,:) = NaN;
    end
end


restind_exone = zeros(cnum-regnum+1,1);
k=2; h=1; 
for n=1:cnum
    if strcmp(cnames{n},list{k})
        if k==regnum
        else
        k=k+1;
        end
    else
        restind_exone(h) = n;
        if h==(cnum-regnum+1)
        else
            h=h+1;
        end
    end
end


newtradematrix = newmat(restind_exone,restind_exone);

% % change country names, and number
cnames_tmp = cnames;
cnames_long_tmp = cnames_long;
cnames_tmp(regind(1)) = name;
cnames_long_tmp(regind(1)) = name_long;

newcnames_tmp = cnames_tmp(restind_exone);
newcnames_long_tmp = cnames_long_tmp(restind_exone);

[newcnames idx] = sort(newcnames_tmp);
newcnames_long = newcnames_long_tmp(idx);

%reorder tradematrix
newtradematrix = newtradematrix(idx,idx);

newcnum = length(newcnames);



