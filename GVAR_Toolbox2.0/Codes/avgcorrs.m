
function [avgcorr avgcorr_d avgcorr_VARXres] = avgcorrs(cnum,cnames,maxlag,dv,endoglist,vnames,epsilon)

%**************************************************************************
% PURPOSE: Calculates Average pair-wise cross-section correlations of all
%          variables and associated model's residuals
%--------------------------------------------------------------------------
% OUTPUT: Nothing, it exports results to 'Output\avgcorrs.xls'
%
% INPUT: 
% - cnum: number of countries
% - cnames: cell list of countries' names (short names)
% - maxlag: maximum lag order of endogenous and weakly exogenous variables
% for each country, i.e., lag order of the GVAR
% - dv: structure, contains data of domestic variables for each country
% - endoglist: structure, contains list of endogenous variables of each
% country (short names)
% - vnames: cell, list of domestic variables (short names)
% - epsilon: struct, contains VARX residuals for each country model
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2010
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

% Creating matrices containing each type of variable or of residual in
% order later to compute correlation matrices.

for j=1:length(vnames)
    
    block.(vnames{j}) = [];
    dblock.(vnames{j}) = [];
    VARXresblock.(vnames{j}) = [];
    
end

for n=1:cnum
    
    k=1;
    for j=1:length(vnames)
        flag = 0;
        for i=1:length(endoglist.(cnames{n}))
            if strcmp(vnames(j),endoglist.(cnames{n})(i))
                flag = 1;
            end
        end
        if flag ==1
            varind.(vnames{j})(n) = 1;
            block.(vnames{j}) = [block.(vnames{j}) dv.(vnames{j}).(cnames{n})(1+maxlag:end)]; % trim series according to estimation sample
            dblock.(vnames{j}) = [dblock.(vnames{j}) diff(dv.(vnames{j}).(cnames{n})(1+maxlag:end))]; % trim series according to estimation sample
            VARXresblock.(vnames{j}) = [VARXresblock.(vnames{j}) epsilon.(cnames{n})(k,:)'];
            k=k+1;
        else
            varind.(vnames{j})(n) = 0;
        end
    end
end

 
% Calculating correlation

for j=1:length(vnames)
    corrmtrx.(vnames{j}) = corrmat(block.(vnames{j}));
    dcorrmtrx.(vnames{j}) = corrmat(dblock.(vnames{j}));
    VARXrescorrmtrx.(vnames{j}) = corrmat(VARXresblock.(vnames{j}));
end

 
% Initialize arrays which will contain average correlations

for j=1:length(vnames)
    avgcorr.(vnames{j}) = NaN(cnum,1);
    avgcorr_d.(vnames{j}) = NaN(cnum,1);
    avgcorr_VARXres.(vnames{j}) = NaN(cnum,1);
end

% Making averages of pairwise correlations

for j=1:length(vnames)
    pointer = 1;
    for n=1:cnum
        if varind.(vnames{j})(n) == 1
            avgcorr.(vnames{j})(n) = (sum(corrmtrx.(vnames{j})(:,pointer))-1)/(rows(corrmtrx.(vnames{j}))-1);
            avgcorr_d.(vnames{j})(n) = (sum(dcorrmtrx.(vnames{j})(:,pointer))-1)/(rows(dcorrmtrx.(vnames{j}))-1);
            avgcorr_VARXres.(vnames{j})(n) = (sum(VARXrescorrmtrx.(vnames{j})(:,pointer))-1)/(rows(VARXrescorrmtrx.(vnames{j}))-1);
            pointer = pointer + 1;
        end
    end
end