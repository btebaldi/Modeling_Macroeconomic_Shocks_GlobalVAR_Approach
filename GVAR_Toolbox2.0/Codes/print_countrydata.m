function print_countrydata(date,cnum,cnames,cnames_long,...
         endoglist,endog,exoglist,exog,gvnames,outdir)

%**************************************************************************
% PURPOSE: printing endogenous and foreign variables for each country
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************

vlab = date;
    for n=1:cnum

        hlab = {'date'};
        tab = [];

        for i=1:length(endoglist.(cnames{n}))
            hlab = [hlab endoglist.(cnames{n})(i)]; %#ok
            tab = [tab endog.(cnames{n})(:,i)]; %#ok
        end

        for j=1:length(exoglist.(cnames{n}))
            fvname_tmp = exoglist.(cnames{n}){j};

            gnameflag = 0;
            if not(isempty(gvnames))
                for g=1:length(gvnames)
                    if strcmp(fvname_tmp,gvnames(g))
                        gnameflag = 1;
                    end
                end
            end

            if gnameflag == 0

                fvname_c = sprintf('%ss',fvname_tmp);
                fvname_s = {fvname_c};
                fvname_tmp = fvname_s;
            end


            hlab = [hlab fvname_tmp]; %#ok
            tab = [tab exog.(cnames{n})(:,j)]; %#ok
        end
        
        toprint = [hlab; vlab num2cell(tab)];

        xlswrite([outdir 'countrydata.xlsx'],toprint,cnames_long{n});        

    end

    