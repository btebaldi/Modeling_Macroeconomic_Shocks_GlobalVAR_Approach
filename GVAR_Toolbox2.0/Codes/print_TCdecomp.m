function print_TCdecomp(x, xp, xc, x_tilde, vnames,gvnames,xnames, date, ...
                        cnames_x,maxlag, outdir)


%**************************************************************************
% PURPOSE: printing TC decomposition output to excel files
%--------------------------------------------------------------------------
%  Meiling He, September 2011
%  CFAP, Judge Business School, University of Cambridge
%  mh590@cam.ac.uk
%**************************************************************************

vlab = date;

vnlist = [vnames gvnames];% create variable list that includes global 
                           % variables and domestic variables.
k = length(xnames);
x_selected_indicator = false(k,length(vnlist));

%write x

for j = 1:length(vnlist)
 
    for i = 1:length(xnames)
        x_selected_indicator(i,j) = strcmp(vnlist(j),xnames(i));
    end
end 

for j = 1:length(vnlist)
 
  
    x_selected = (x(x_selected_indicator(:,j),:))';
    
    hlab = ['Date' (cnames_x(x_selected_indicator(:,j)))'];
    
    toprint = [hlab; vlab num2cell(x_selected)];
    sheet_name = vnlist{j};
    
    xlswrite([outdir 'TCdecomposition.xls'],toprint,sheet_name);
    
end 

%write xp 
for j = 1:length(vnlist)
 
       
    xp_selected = (xp(x_selected_indicator(:,j),:))';
    
    hlab = ['Date' (cnames_x(x_selected_indicator(:,j)))'];
    nc = sum(x_selected_indicator(:,j));
    if maxlag ~= 0
        xp_selected = [NaN(maxlag,nc); xp_selected]; %#ok
    end
    
    toprint = [hlab; vlab num2cell(xp_selected)];
    sheet_name = [vnlist{j} '_p'];
    
    xlswrite([outdir 'TCdecomposition.xls'],toprint,sheet_name);
    
end 

%write xc 
for j = 1:length(vnlist)
 
       
    xc_selected = (xc(x_selected_indicator(:,j),:))';
    
    hlab = ['Date' (cnames_x(x_selected_indicator(:,j)))'];
    nc = sum(x_selected_indicator(:,j));
    if maxlag ~= 0
        xc_selected = [NaN(maxlag,nc); xc_selected]; %#ok
    end
    
    toprint = [hlab; vlab num2cell(xc_selected)];
    sheet_name = [vnlist{j} '_c'];
    
    xlswrite([outdir 'TCdecomposition.xls'],toprint,sheet_name);
    
end 

% for j = 1:length(vnlist)
%  
%        
%     x_tilde_selected = (x_tilde(x_selected_indicator(:,j),:))';
%     
%     hlab = ['Date' (cnames_x(x_selected_indicator(:,j)))'];
%     nc = sum(x_selected_indicator(:,j));
%     if maxlag ~= 0
%         x_tilde_selected = [NaN(maxlag,nc); x_tilde_selected]; %#ok
%     end
%     
%     toprint = [hlab; vlab num2cell(x_tilde_selected)];
%     sheet_name = [vnlist{j} '_dev'];
%     
%     xlswrite([outdir 'TCdecomposition.xls'],toprint,sheet_name);
%     
% end 
        
    
    
   