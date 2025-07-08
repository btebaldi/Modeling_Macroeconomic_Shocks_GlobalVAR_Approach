
function fevddir = print_fevds(cnames,cnames_long,rnames,rnames_long,xnames,cnames_x,N,FEVD, ...
                    vnames,vnames_long,gvnames,gvnames_long,shocktype,shocksign,...
                    vartype,countype,regtype,outdir_t,bsflag,regflag,sgirfflag,label1)

%**************************************************************************
% PURPOSE: printing generalised forecast error variance decomposition
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************


if sgirfflag == 0
    title = {'Generalized Forecast Error Variance Decompositions'};
    foldername = {'GFEVDs'};
elseif sgirfflag == 1
    title = {'Structural Generalized Forecast Error Variance Decompositions'};
    foldername = {'SGFEVDs'};
elseif sgirfflag == 2
    title = {'Orthogonalized Forecast Error Variance Decompositions'};
    foldername = {'OFEVDs'};
end


allvlist = [vnames gvnames];
allvlist_long = [vnames_long gvnames_long];

for i=1:length(allvlist)
    if strcmp(allvlist{i},vartype)
        vartype_long = allvlist_long(i);
    end
end

for n=1:length(cnames)
    if strcmp(cnames{n},countype)
        countype_long = cnames_long(n);
    end
end

for r=1:length(rnames)
    if strcmp(rnames{r}, regtype)
        regtype_long = rnames_long(r);
    end
end



if shocktype == 3
    if shocksign == -1
        fevddir = sprintf('Output\\%s\\%s\\Global 1se neg shock to %s\\',outdir_t,foldername{1},vartype_long{1});
    elseif shocksign == 1
        fevddir = sprintf('Output\\%s\\%s\\Global 1se pos shock to %s\\',outdir_t,foldername{1},vartype_long{1});
    end
elseif shocktype == 2
    if shocksign == -1
        fevddir = sprintf('Output\\%s\\%s\\%s 1se neg shock to %s\\',outdir_t,foldername{1},regtype_long{1},vartype_long{1});
    elseif shocksign == 1
        fevddir = sprintf('Output\\%s\\%s\\%s 1se pos shock to %s\\',outdir_t,foldername{1},regtype_long{1},vartype_long{1});
    end
elseif shocktype == 1
    if shocksign == -1
        fevddir = sprintf('Output\\%s\\%s\\%s 1se neg shock to %s\\',outdir_t,foldername{1},countype_long{1},vartype_long{1});
    elseif shocksign == 1
        fevddir = sprintf('Output\\%s\\%s\\%s 1se pos shock to %s\\',outdir_t,foldername{1},countype_long{1},vartype_long{1});
    end
end

mkdir(fevddir);






if shocktype == 3 % global shock
    subtitle = sprintf('Proportion of the N-step ahead Forecast Error Variance of %s Explained by Conditioning on Contemporaneous and Future Innovations of the Country Equations', vartype_long{1});
elseif shocktype == 2 % rs shock
    subtitle = sprintf('Proportion of the N-step ahead Forecast Error Variance of %s %s Explained by Conditioning on Contemporaneous and Future Innovations of the Country Equations', regtype_long{1},vartype_long{1});
elseif shocktype == 1 % cs shock
    subtitle = sprintf('Proportion of the N-step ahead Forecast Error Variance of %s %s Explained by Conditioning on Contemporaneous and Future Innovations of the Country Equations', countype_long{1},vartype_long{1});
end

subtitle = {subtitle};

timetrend = 0:N;



if regflag == 0
    obj1 = [[title {''}; subtitle {''}; {''} {''}];{''} {''}; [cnames_x xnames]];
    obj2 = [timetrend; FEVD];
else
    rlist_tot=[];
    labelvar_tot=[];
    totalblock = [];
    for i=1:length(allvlist)
        if isfield(FEVD,allvlist{i})
            rlist = fieldnames(FEVD.(allvlist{i}));
            
            rFblock = [];
            labelvar = [];
            for j=1:length(rlist)
                rFblock = [rFblock; FEVD.(allvlist{i}).(rlist{j})]; %#ok
                labelvar = [labelvar; allvlist(i)]; %#ok
            end
            
            rlist_tot = [rlist_tot; rlist]; %#ok
            labelvar_tot = [labelvar_tot; labelvar]; %#ok
            
            varblock = rFblock;
            totalblock = [totalblock; varblock]; %#ok
        end
    end
    
    obj1 = [[title {''}; subtitle {''}; {''} {''}];{''} {''}; [rlist_tot labelvar_tot]];
    obj2 = [timetrend; totalblock];
end

warning off MATLAB:xlswrite:AddSheet




if regflag == 0
    if bsflag == 0
        xlswrite('Output\fevds.xls',obj1,label1,'A1');
        xlswrite('Output\fevds.xls',obj2,label1,'C4');
    else
        xlswrite('Output\fevds_bs.xls',obj1,label1,'A1');
        xlswrite('Output\fevds_bs.xls',obj2,label1,'C4');
    end
else
    if bsflag == 0
        xlswrite('Output\rfevds.xls',obj1,label1,'A1');
        xlswrite('Output\rfevds.xls',obj2,label1,'C4');
    else
        xlswrite('Output\rfevds_bs.xls',obj1,label1,'A1');
        xlswrite('Output\rfevds_bs.xls',obj2,label1,'C4');
    end
end

