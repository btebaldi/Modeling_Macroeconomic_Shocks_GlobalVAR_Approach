function print_VARX_SC(cnum,cnames,cnames_long,aic_sbc,F_serialcorr,vnames,gvnames,vnum,gvnum,dvflag,gvflag,outdir)


%**************************************************************************
% PURPOSE: printing VARX* statistics
%--------------------------------------------------------------------------
% Alessandro Galesi, July 2011
% Centro de Estudios Monetarios y Financieros, Madrid.
% alessandro.galesi@cemfi.edu.es
%**************************************************************************     
     
tab = ['Choice Criteria for Selecting the Order of the VARX* Models Together With Corresponding Residual Serial Correlation F-Statistics' num2cell(NaN(1,8+vnum+gvnum))];
tab = [tab; num2cell(NaN(1,9+vnum+gvnum))];
tab = [tab; {' ' 'p' 'q' 'AIC' 'SBC' 'logLik' ' ' ' ' 'Fcrit_0.05'} vnames gvnames];


for n=1:cnum
    for j=1:rows(aic_sbc.(cnames{n}))
        clab = cnames_long(n); 
        out_tmp = aic_sbc.(cnames{n})(j,:);
        out1 = [out_tmp(4:5) out_tmp(2:3) out_tmp(1)]; % reorder results, only esthetical purpose
        
        flab_tmp = sprintf('F(%d,%d)', F_serialcorr.(cnames{n})(j,1), F_serialcorr.(cnames{n})(j,2));
        flab = {flab_tmp};
        vec = flab; %#ok
        
        critfigs = F_serialcorr.(cnames{n})(j,3);
        dvfigs = zeros(1,vnum+gvnum);
        k=1;
        for i=1:vnum
            if dvflag(n,i) == 1
                dvfigs(i) = F_serialcorr.(cnames{n})(j,k+3);
                k=k+1;
            else
                dvfigs(i) = NaN;
            end
        end
        
        if not(gvnum == 0)
            for i=1:gvnum
                if gvflag(n,i) == 2
                    dvfigs(length(vnames)+i) = F_serialcorr.(cnames{n})(j,k+3);
                    k=k+1;
                else
                    dvfigs(length(vnames)+i) = NaN;
                end
            end
        end
        out2 = [critfigs dvfigs];
        
        tab = [tab; clab num2cell(out1) ' ' flab num2cell(out2)]; %#ok
    end
end

xlswrite([outdir 'output.xlsx'],tab,'VARX_sc');


