function plot_irfs(vnames,gvnames,gxvnames,vnames_long,gvnames_long,vartype,rnames,regtype,rnames_long,cnames,countype,cnames_long,...
         irfdir,shocktype,shocksign,N,K,cnames_x,xnames,IRF,label1,label2,copyfile_flag,rgraphs_flag,sgirfflag)

%**************************************************************************
% PURPOSE: plotting impulse responses 
%--------------------------------------------------------------------------
% Alessandro Galesi, 2014. CEMFI, Madrid. galesi@cemfi.edu.es
%**************************************************************************    

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


for i=1:length(allvlist_long)
    
   
    vlabel_tmp = allvlist_long{i};
    
    if rgraphs_flag == 0
        xlsoutfile = sprintf('%s %s.xls',label1,vlabel_tmp);
    elseif rgraphs_flag == 1
        xlsoutfile = sprintf('r%s %s.xls',label1,vlabel_tmp);
        
        % avoid global exogenous variables (which belong to the dominant
        % unit model)
        flag = 0;
        for j=1:length(gxvnames);
            if strcmp(allvlist{i},gxvnames{j})
                flag = 1;
            end
        end        
        
    end
    xlsoutdir = [irfdir xlsoutfile];

    if rgraphs_flag == 1 && flag == 1
        continue
    end
    
    
    if copyfile_flag == 1 % copy the template file
        message = sprintf('creating %s',xlsoutfile);
        disp(message);         
        
        templatefile = sprintf('Tech\\irfs_%s.xls',label1);
        copyfile(templatefile,xlsoutdir);
    end

    
    if sgirfflag == 0
        title = {'Generalized Impulse Response Functions'};
    elseif sgirfflag == 1
        title = {'Structural Generalized Impulse Response Functions'};
    elseif sgirfflag == 2
        title = {'Orthogonalized Impulse Response Functions'};
    end

    if shocktype == 3
        if shocksign == 1
            subtitle = sprintf('One s.e. Positive Global Shock to %-1s', vartype_long{1});
        elseif shocksign == -1
            subtitle = sprintf('One s.e. Negative Global Shock to %-1s', vartype_long{1});
        end
    elseif shocktype == 2
        if shocksign == 1
            subtitle = sprintf('One s.e. Positive Region-specific shock to %s %s', regtype_long{1},vartype_long{1});
        elseif shocksign == -1
            subtitle = sprintf('One s.e. Negative Global Shock to %s %s', regtype_long{1}, vartype_long{1});
        end
    elseif shocktype == 1
        if shocksign == 1
            subtitle = sprintf('One s.e. Positive Shock to %s %s', countype_long{1},vartype_long{1});
        elseif shocksign == -1
            subtitle = sprintf('One s.e. Negative Shock to %s %s', countype_long{1},vartype_long{1});
        end
    end

    subtitle = {subtitle};

    subsubtitle1 = label2;

    timetrend = 0:N;


    % retrieve variable-specific IRFs
    
    if rgraphs_flag == 0  % country-specific irfs
        
    IRF_sp = [];
    cnames_x_sp = {};
    xnames_sp = {};

    for s=1:K
        if strcmp(xnames(s),allvlist(i))
            IRF_sp = [IRF_sp; IRF(s,:)]; %#ok
            cnames_x_sp = [cnames_x_sp; cnames_x(s)]; %#ok
            xnames_sp = [xnames_sp; xnames(s)]; %#ok
        end
    end
    
    obj1 = [[title {''}; subtitle {''}; subsubtitle1 {''}];{''} {''}; [cnames_x_sp xnames_sp]];
    obj2 = [timetrend; IRF_sp];
    
    elseif rgraphs_flag == 1 % regional irfs
        
        if isfield(IRF,allvlist{i})
            rlist = fieldnames(IRF.(allvlist{i}));
            
            rIRF_sp = [];
            labelvar = [];
            for j=1:length(rlist)
                rIRF_sp = [rIRF_sp; IRF.(allvlist{i}).(rlist{j})]; %#ok
                labelvar = [labelvar; allvlist(i)]; %#ok
            end
            
            obj1 = [[title {''}; subtitle {''}; subsubtitle1 {''}];{''} {''}; [rlist labelvar]];
            obj2 = [timetrend; rIRF_sp];
        end
    end
    
    warning off MATLAB:xlswrite:AddSheet
   
    % print point estimates
    xlswrite(xlsoutdir,obj1,label2,'A1');
    xlswrite(xlsoutdir,obj2,label2,'C4');
    
end
