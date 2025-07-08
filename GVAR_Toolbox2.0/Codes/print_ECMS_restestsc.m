function print_ECMS_restestsc(cnum,cnames,cnames_long,psc_degfrsc_Fcrit_Fsc,vnames,gvnames,vnum,gvnum,dvflag,gvflag,outdir)

%**************************************************************************
% PURPOSE: printing serial correlation test statistics of VECMX* residuals
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
    

title = {'F-Statistics for the Serial Correlation Test of the VECMX* Residuals'};

clab = [];
for n=1:length(cnames_long)
    clab = [clab; cnames_long(n)]; %#ok
end

% print FSC serial correlation

vec = [];
for n = 1:cnum
    for j=1:rows(psc_degfrsc_Fcrit_Fsc.(cnames{n}))
        flab_tmp = sprintf('F(%d,%d)', psc_degfrsc_Fcrit_Fsc.(cnames{n})(j,1), psc_degfrsc_Fcrit_Fsc.(cnames{n})(j,2));
        flab = {flab_tmp};
        vec = [vec; flab]; %#ok
    end
end

hlab = ['Fcrit_0.05' vnames gvnames];

table = [];
for n = 1:cnum
    for j=1:rows(psc_degfrsc_Fcrit_Fsc.(cnames{n}))
        
        critfigs = psc_degfrsc_Fcrit_Fsc.(cnames{n})(j,3);
        dvfigs = zeros(1,vnum+gvnum);
        k=1;
        for i=1:vnum
            if dvflag(n,i) == 1
                dvfigs(i) = psc_degfrsc_Fcrit_Fsc.(cnames{n})(j,k+3);
                k=k+1;
            else
                dvfigs(i) = NaN;
            end
        end
        
        if not(gvnum == 0)
            for i=1:gvnum
                if gvflag(n,i) == 2
                    dvfigs(length(vnames)+i) = psc_degfrsc_Fcrit_Fsc.(cnames{n})(j,k+3);
                    k=k+1;
                else
                    dvfigs(length(vnames)+i) = NaN;
                end
            end
        end
        
        table = [table; critfigs dvfigs]; %#ok
    end
end


tab = [title num2cell(NaN(1,1+cols(table)))];
tab = [tab; num2cell(NaN(1,2+cols(table)))];
tab = [tab; num2cell(NaN(1,2)) hlab];
tab = [tab; clab vec num2cell(table)];    

xlswrite([outdir 'output.xlsx'],tab,'ECMS_restestsc');

