function print_forecasts(lastobs,fhorz,cnames_s_x,xnames,xforc,rx,outdir)

%**************************************************************************
% PURPOSE: printing forecasts of GVAR
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     

ldate = lastobs;
if not(iscell(lastobs))
    ldate_tmp = num2str(ldate);
    ldate = {ldate_tmp};
    % annual data
end
stitle1 = ['Last in-sample estimation date:' ldate];
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

stitle2 = ['Last forecast date:' fdate(end)];

% create labels
vtitle = {};
for i=1:rows(xforc)
    vtitlea = [cnames_s_x(i) xnames(i) 'Forecast'];
    vtitleb = {' ' ' ' 'Actual'};
    vtitle = [vtitle; vtitlea; vtitleb]; %#ok
end

% create table of data
obsdiff = abs(cols(xforc)-cols(rx));

nanbar = [];
for i=1:obsdiff
    nanbar = [nanbar NaN]; %#ok
end

table = [];
for i=1:rows(xforc)
    if cols(xforc)>cols(rx) && not(isempty(rx))
        table = [table; xforc(i,:); rx(i,:) nanbar]; %#ok
    elseif cols(xforc)>cols(rx) && isempty(rx)
        table = [table; xforc(i,:); nanbar]; %#ok
    else
        table = [table; xforc(i,:) nanbar; rx(i,:)]; %#ok
    end
end


tab = ['GVAR Forecasts' num2cell(NaN(1,2+cols(table)))];
tab = [tab; stitle1 num2cell(NaN(1,1+cols(table)))];
tab = [tab; stitle2 num2cell(NaN(1,1+cols(table)))];
tab = [tab; num2cell(NaN(1,3+cols(table)))];
tab = [tab; num2cell(NaN(1,3)) fdate'];
tab = [tab; vtitle num2cell(table)];

xlswrite([outdir 'output.xlsx'],tab,'forecasts');

