function [TC_Restr] = TC_trend_restr(cnames_x, xnames,outdir) 

%**************************************************************************
% PURPOSE: It allows the user to input trend restrictions for TC decomposition. 
%--------------------------------------------------------------------------
% Meiling He and L. Vanessa Smith, September 2011
% CFAP, Judge Business School, University of Cambridge
% mh590@cam.ac.uk and lvs21@cam.ac.uk
%**************************************************************************
k = length(cnames_x);

title = {['Insert Trend Restrictions for the Trend/Cycle Decomposition of the GVAR '] [''] ['']}; %#ok
hlab = {['Country'] ['Variable'] ['Trend Restrictions']}; %#ok
hlab = [cell(1, cols(hlab)) ; hlab ];
empty_cell = cell(k,1);

tab = [title; hlab;cnames_x xnames empty_cell];

xlswrite([outdir 'TC_trend_restr.xls'],tab,'TC_trend_restr');

% pause and let user input trend restrictions
%**********************************************************************
winopen([outdir 'TC_trend_restr.xls']);
disp(' ');
disp('>>> Pause and go to TC_trend_restr.xls: Insert trend restictions for the Trend/Cycle(TC) decomposition,');
disp('    then press enter.');
disp(' ');
disp('======================================================================================');
disp('Note:');
disp('Insert the value of 1 next to the country variable that you wish to impose a trend');
disp('restriction on, leaving all other cells empty.');

disp(' ');
disp(' USING THE FULL DEMO INTERFACE FILE');
disp('If you are using the full demo interface file to perform the TC decomposition with');
disp('trend restrictions, restrict the trend for inflation and the interest rate (short'); 
disp('and long) for all countries.');
disp('=======================================================================================');
disp(' ');
pause
disp('>>> Make sure you have saved and closed the TC_trend_restr.xls file. If so, press enter again.');
endpause = input(' ','s'); %#ok
disp('- The program is now running (do not press any key)');
disp(' ');

% read inputted trend restrictions
%********************************************

TC_Restr = false(k,1);

TC_Restr_temp  = xlsread([outdir 'TC_trend_restr.xls'],'TC_trend_restr',['C4:C' num2str(3+k)]);

if length(TC_Restr_temp)~=k && isempty(TC_Restr_temp)~=1 && sum(TC_Restr_temp)~=0
    
    for start_ind = 1: k

        temp  = xlsread([outdir 'TC_trend_restr.xls'],'TC_trend_restr',['C' num2str(3+start_ind)]);

        if ~isempty(temp)
            break;
        end
    end 

    end_ind = k;
    while end_ind >start_ind 

        temp  = xlsread([outdir 'TC_trend_restr.xls'],'TC_trend_restr',['C' num2str(3+end_ind)]);

        if ~isempty(temp)
            break;
        end
        end_ind=end_ind-1;
    end

    if start_ind<=end_ind

        TC_Restr_temp  = xlsread([outdir 'TC_trend_restr.xls'],'TC_trend_restr',['C' num2str(3+start_ind) ':C' num2str(3+end_ind)]);

        TC_Restr(start_ind:end_ind,1) = ~isnan(TC_Restr_temp)&(TC_Restr_temp~=0);

    end

elseif length(TC_Restr_temp)==k && isempty(TC_Restr_temp)~=1 && sum(TC_Restr_temp)~=0
    TC_Restr = ~isnan(TC_Restr_temp)&(TC_Restr_temp~=0);
end 


