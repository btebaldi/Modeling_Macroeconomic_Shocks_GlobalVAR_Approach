
function irfdir = print_irfs(outdir_t,cnames,cnames_long,rnames,rnames_long,vnames,gvnames,vnames_long,gvnames_long,xnames,cnames_x,N,IRF, ...
                    shocktype,shocksign,vartype,countype,regtype,bsflag,regflag,sgirfflag,label1)

%**************************************************************************
% PURPOSE: prints impulse response functions
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************   
     

if sgirfflag == 0
    title = {'Generalized Impulse Response Functions'};
    foldername = {'GIRFs'};
elseif sgirfflag == 1
    title = {'Structural Generalized Impulse Response Functions'};
    foldername = {'SGIRFs'};
elseif sgirfflag == 2
    title = {'Orthogonalized Impulse Response Functions'};
    foldername = {'OIRFs'};
end

allvlist = [vnames gvnames];
allvlist_long = [vnames_long gvnames_long];

for i=1:length(allvlist)
    if strcmp(allvlist{i},vartype)
        vartype_long = allvlist_long(i);
    end
end

for r=1:length(rnames)
    if strcmp(rnames{r}, regtype)
        regtype_long = rnames_long(r);
    end
end

for n=1:length(cnames)
    if strcmp(cnames{n},countype)
        countype_long = cnames_long(n);
    end
end



if shocktype == 3
    if shocksign == -1
        irfdir = sprintf('Output\\%s\\%s\\Global 1se neg shock to %s\\',outdir_t,foldername{1},vartype_long{1});
    elseif shocksign == 1
        irfdir = sprintf('Output\\%s\\%s\\Global 1se pos shock to %s\\',outdir_t,foldername{1},vartype_long{1});
    end
elseif shocktype == 2
    if shocksign == -1
        irfdir = sprintf('Output\\%s\\%s\\%s 1se neg shock to %s\\',outdir_t,foldername{1},regtype_long{1},vartype_long{1});
    elseif shocksign == 1
        irfdir = sprintf('Output\\%s\\%s\\%s 1se pos shock to %s\\',outdir_t,foldername{1},regtype_long{1},vartype_long{1});
    end
elseif shocktype == 1
    if shocksign == -1
        irfdir = sprintf('Output\\%s\\%s\\%s 1se neg shock to %s\\',outdir_t,foldername{1},countype_long{1},vartype_long{1});
    elseif shocksign == 1
        irfdir = sprintf('Output\\%s\\%s\\%s 1se pos shock to %s\\',outdir_t,foldername{1},countype_long{1},vartype_long{1});
    end
end

mkdir(irfdir);




if shocktype == 3
    if shocksign == 1
        subtitle = sprintf('One s.e. Positive Global Shock to %-1s', vartype_long{1});
    elseif shocksign == -1
        subtitle = sprintf('One s.e. Negative Global Shock to %-1s', vartype_long{1});
    end
elseif shocktype == 2
    if shocksign == 1
        subtitle = sprintf('One s.e. Positive Regional Shock to %s %s', regtype_long{1},vartype_long{1});
    elseif shocksign == -1
        subtitle = sprintf('One s.e. Negative Regional Shock to %s %s', regtype_long{1},vartype_long{1});
    end
elseif shocktype == 1
    if shocksign == 1
        subtitle = sprintf('One s.e. Positive Shock to %s %s', countype_long{1},vartype_long{1});
    elseif shocksign == -1
        subtitle = sprintf('One s.e. Negative Shock to %s %s', countype_long{1},vartype_long{1});
    end
end

subtitle = {subtitle};

subsubtitle1 = {label1};

timetrend = 0:N;


if regflag == 0
    obj1 = [[title {''}; subtitle {''}; subsubtitle1 {''}];{''} {''}; [cnames_x xnames]];
    obj2 = [timetrend; IRF];
else
    rlist_tot=[];
    labelvar_tot=[];
    totalblock = [];
    for i=1:length(allvlist)
        if isfield(IRF,allvlist{i})
            rlist = fieldnames(IRF.(allvlist{i}));
            
            rIblock = [];
            labelvar = [];
            for j=1:length(rlist)
                rIblock = [rIblock; IRF.(allvlist{i}).(rlist{j})]; %#ok
                labelvar = [labelvar; allvlist(i)]; %#ok
            end
            
            rlist_tot = [rlist_tot; rlist]; %#ok
            labelvar_tot = [labelvar_tot; labelvar]; %#ok
            
            varblock = rIblock;
            totalblock = [totalblock; varblock]; %#ok
        end
    end
    
    obj1 = [[title {''}; subtitle {''}; subsubtitle1 {''}];{''} {''}; [rlist_tot labelvar_tot]];
    obj2 = [timetrend; totalblock];
    
end

warning off MATLAB:xlswrite:AddSheet


if regflag == 0
    if bsflag == 0
        xlswrite('Output\irfs.xls',obj1,label1,'A1');
        xlswrite('Output\irfs.xls',obj2,label1,'C4');
    else
        xlswrite('Output\irfs_bs.xls',obj1,label1,'A1');
        xlswrite('Output\irfs_bs.xls',obj2,label1,'C4');
    end
else
    if bsflag == 0
        xlswrite('Output\rirfs.xls',obj1,label1,'A1');
        xlswrite('Output\rirfs.xls',obj2,label1,'C4');
    else
        xlswrite('Output\rirfs_bs.xls',obj1,label1,'A1');
        xlswrite('Output\rirfs_bs.xls',obj2,label1,'C4');
    end
end
