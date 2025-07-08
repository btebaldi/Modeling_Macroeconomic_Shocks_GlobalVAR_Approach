function [lb_restr] = ExF_lb_restr(cnames_x, xnames,outdir) 

%**************************************************************************
% PURPOSE: It allows the user to input lower bound restrictions for the
%          computation of GVAR ex-ante forecasts. 
%--------------------------------------------------------------------------
% L. Vanessa Smith, Cambridge, January 2014
% lvs21@cam.ac.uk
%**************************************************************************


title = {['Insert Lower Bound Restrictions for the Ex-ante Forecasts of the GVAR '] [''] ['']}; %#ok
hlab = {['Country'] ['Variable'] ['Lower Bound Restrictions']}; %#ok
hlab = [cell(1, cols(hlab)) ; hlab ];
empty_cell = cell(length(cnames_x),1);

tab = [title; hlab;cnames_x xnames empty_cell];

xlswrite([outdir 'ExF_lb_restr.xls'],tab,'restrictions');


% pause and let user input trend restrictions
%**********************************************************************
winopen([outdir 'ExF_lb_restr.xls']);
disp(' ');
disp('>>> Pause and go to ExF_lb_restr.xls: Insert lower bound restrictions for the country-specific');
disp('    variables of interest, then press enter.');
disp(' ');
pause
disp('>>> Make sure you have saved and closed the ExF_lb_restr.xls file. If so, press enter again.');
endpause = input(' ','s'); %#ok
disp('- The program is now running (do not press any key)');
disp(' ');

% read inputted lower bound restrictions
%********************************************

k=length(cnames_x);

lb_restr = NaN(k,1);


lb_restr_temp  = xlsread([outdir 'ExF_lb_restr.xls'],'restrictions',['C4:C' num2str(3+k)]);

if length(lb_restr_temp)~=k 
    
    for start_ind = 1: k

        temp  = xlsread([outdir 'ExF_lb_restr.xls'],'restrictions',['C' num2str(3+start_ind)]);

        if ~isempty(temp)
            break;
        end
    end 

    end_ind = k;
    while end_ind >start_ind 

        temp  = xlsread([outdir 'ExF_lb_restr.xls'],'restrictions',['C' num2str(3+end_ind)]);

        if ~isempty(temp)
            break;
        end
        end_ind=end_ind-1;
    end

    if start_ind<=end_ind

        lb_restr_temp  = xlsread([outdir 'ExF_lb_restr.xls'],'restrictions',['C' num2str(3+start_ind) ':C' num2str(3+end_ind)]);
 
        lb_restr(start_ind:end_ind,1) = lb_restr_temp;

    end

elseif length(lb_restr_temp)==k 

 lb_restr = lb_restr_temp;
 
end 
