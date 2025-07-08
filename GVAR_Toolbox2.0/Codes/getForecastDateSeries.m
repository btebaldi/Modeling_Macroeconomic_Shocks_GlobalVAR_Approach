function [fdate]= getForecastDateSeries(lastobs,fhorz)

%**************************************************************************
% PURPOSE: creat the Dates for the forecast horizons
%--------------------------------------------------------------------------
% Meiling He, Nov 2012
% Faculty of Economics, University of Cambridge
% mh590@cam.ac.uk
%**************************************************************************
%

ldate = lastobs;
if not(iscell(lastobs))
    ldate_tmp = num2str(ldate);
    ldate = {ldate_tmp};
    % annual data
end
% stitle1 = ['Last in-sample estimation date:' ldate];
ldate_ch = char(ldate);

ysw = [];
if not(iscell(lastobs)); % annual data
    freq = 1;
    
elseif strcmp(ldate_ch(5),'Q'); % quarterly data
    freq = 4;
    for i=1:freq
        qtmp = sprintf('%d',i);
        if strcmp(ldate_ch(6),qtmp)
            if i==freq
                foreq = 1;
                ysw = 1;
            else
                foreq = i+1;
                ysw = 0;
            end
        end
    end
elseif strcmp(ldate_ch(5),'M'); % monthly data
    freq = 12;
    for i=1:freq
        if strcmp(ldate_ch(6),'0')
            qtmp = sprintf('0%d',i);
        else
            qtmp = sprintf('%d',i);
        end
        if strcmp(ldate_ch(6:7),qtmp)
            if i==freq
                forem = 1;
                ysw = 1;
            else
                forem = i+1;
                ysw = 0;
            end
        end
    end
end

% read the year
if freq == 1 % annual
else % quarterly and monthly
    yr = cellstr(ldate_ch(1:4));
    ymin=1000;
    ymax=2200;
    for i=ymin:ymax
        ytmp = sprintf('%d',i);
        if strcmp(yr,ytmp)
            forey = i+ysw;
        end
    end
end

% create forecast dates
fdate = [];
if freq == 1
    for i=1:fhorz
        fd_num = lastobs+i;
        fd_tmp = num2str(fd_num);
        fd = {fd_tmp};
        fdate = [fdate; fd]; %#ok
    end
elseif freq == 4
    fmin = sprintf('%dQ%d',forey,foreq);
    fdate = {fmin};
    for i=2:fhorz
        if foreq == freq
            forey = forey+1;
            foreq = 1;
        else
            foreq = foreq+1;
        end
        fd = sprintf('%dQ%d',forey,foreq);
        fdate = [fdate; fd]; %#ok
    end
elseif freq == 12
    fmin = sprintf('%dM%d',forey,forem);
    fdate = {fmin};
    for i=2:fhorz
        if forem == freq
            forey = forey+1;
            forem = 1;
        else
            forem = forem+1;
        end
        fd = sprintf('%dM%d',forey,forem);
        fdate = [fdate; fd]; %#ok
    end
end

% stitle2 = ['Last forecast date:' fdate(end)];