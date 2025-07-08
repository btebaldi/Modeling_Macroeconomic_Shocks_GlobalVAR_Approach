function [con_forc_restr_mtx,con_forc_restr_x] = con_forc_restr(K,con_fhorz_restr,xnames,lastobs,cnames_s_x,outdir)

%**************************************************************************
% PURPOSE: creat the xls file for conditional forecast
%--------------------------------------------------------------------------
% Meiling He, Nov 2012
% Faculty of Economics, University of Cambridge
% mh590@cam.ac.uk
%**************************************************************************


%Get date series for forecast
con_fdate = getForecastDateSeries(lastobs,con_fhorz_restr);

horz_lab = [{'Country'},{'Variable'},con_fdate'];

vert_lab = [cnames_s_x,xnames];

empty_mtrx = num2cell(NaN(rows(vert_lab),con_fhorz_restr));

tab = [horz_lab;vert_lab,empty_mtrx];

xlswrite([outdir 'con_forc_restr.xls'],tab,'restrictions');

winopen([outdir 'con_forc_restr.xls']);

disp(' ');
disp('>>> Pause and go to con_forc_restr.xls: Input the conditional forecast restrictions on the');
disp('    variables of interest making sure that they are imposed over the entire restriction');
disp('    horizon, then press enter.');
disp(' ');
disp('=========================================================================================');
disp(' USING THE FULL DEMO INTERFACE FILE');
disp('If you are using the full demo interface file, and would like to replicate the results');
disp('in the Output Full Demo folder:                           ');
disp('1. For the US short-term interest rate (r):');
disp('Set the value of 0.01 for each quarter of the restriction horizon 2013Q2-2014Q1.');
disp('2. For the US long-term interest rate (lr):');
disp('Set the value of 0.02 for each quarter of the restriction horizon 2013Q2-2014Q1.');
disp('=========================================================================================');
disp(' ');
pause
disp('>>> Make sure you have saved and closed the con_forc_restr.xls file. If so, press enter again.');
endpause = input(' ','s'); %#ok
disp('- The program is now running (do not press any key)');
disp(' ');


temp_col_begin = 2+1;
temp_col_end = 2+con_fhorz_restr;

temp_row_begin = 1+1;
temp_row_end = 1+K;

Alphabet=char('a'+(1:26)-1)';

[I,J]=ndgrid(1:26,1:26); 
I=I'; J=J'; 

XX=[Alphabet(I(:)), Alphabet(J(:))];
XX=strvcat(Alphabet,XX); %#ok


range_start = [XX(temp_col_begin,1),num2str( temp_row_begin)];

if temp_col_end<=26
    range_end = [XX(temp_col_end,1),num2str(temp_row_end)];
elseif temp_col_end>26
    range_end = [XX(temp_col_end,:),num2str(temp_row_end)];
end

dummy_check = 1;
while dummy_check==1
    
    % - Substitute this line:
%     [con_forc_restr_mtx,junk,con_forc_restr_mtx_raw]= xlsread([outdir 'con_forc_restr.xls'],'restrictions',[range_start,':',range_end]); %#ok         
    % - with the following ones:
    %---------------------------------------------------------
    [con_forc_restr_mtx_tmp,junk,con_forc_restr_mtx_raw]= xlsread([outdir 'con_forc_restr.xls'],'restrictions',[range_start,':',range_end]); %#ok    
    con_forc_restr_mtx = [];
    for j=1:rows(con_forc_restr_mtx_tmp)
        if not(isnan(sum(con_forc_restr_mtx_tmp(j,:))))
            con_forc_restr_mtx = [con_forc_restr_mtx; con_forc_restr_mtx_tmp(j,:)]; %#ok
        end
    end
    clear con_forc_restr_mtx_tmp
    %----------------------------------------------------------
    % - this solves bug which arises when restrictions are imposed on
    % variables which are noncontiguous (in terms of ordering of the
    % variables in the GVAR). A.G. 2013
    
    
    con_forc_restr_mtx_emptyind = cellfun(@isnan,con_forc_restr_mtx_raw);

    sum_emptyind = sum(con_forc_restr_mtx_emptyind,2);

    if sum(sum_emptyind>0&sum_emptyind<con_fhorz_restr)>0
        disp('>>> Please make sure that the restrictions are imposed across all cells up until the end')      
        disp('    of the restriction horizon. ')
        pause
        disp(' ');
        winopen([outdir 'con_forc_restr.xls']);

       disp(' ');
    disp('>>> Pause and go to con_forc_restr.xls: Input the conditional forecast restrictions on the');
    disp('    variables of interest making sure that they are imposed over the entire restriction');
    disp('    horizon, then press enter.');
    disp(' ');
    disp('=========================================================================================');
    disp(' USING THE FULL DEMO INTERFACE FILE');
    disp('If you are using the full demo interface file, and would like to replicate the results');
    disp('in the Output Full Demo folder downloaded with the program:                           ');
    disp('1. For the US short-term interest rate (r):');
    disp('Set the value of 0.01 quarter of the restriction horizon, namely quarters 2013Q2-2014Q1.');
    disp('2. For the US long-term interest rate (lr):');
    disp('Set the value of 0.02 quarter of the restriction horizon, namely quarters 2013Q2-2014Q1.');
    disp('=========================================================================================');
    disp(' ');
    pause
    disp('>>> Make sure you have saved and closed the con_forc_restr.xls file. If so, press enter again.');
    endpause = input(' ','s'); %#ok
    disp('- The program is now running (do not press any key)');
    disp(' ');

    else
        dummy_check=0;
    end 
end 

con_forc_restr_x = (sum_emptyind==0);

    