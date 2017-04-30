%% Subroutine of MAIN.
%--------------Perform Fit--------------------------
% Fitting assigns values to Xsol, Z, deltaR_model, ratio_model, LCTE.
% LCTE changes if eta changes by changing L, or if the entire LCTE
% changes from self-consistent temperature changes.
if manualfit
    if ~doughnut
        switch sigfit
            case 1, fprintf('Manual fitting to V(in)(f) normalized at %0.1f MHz...\n',f_data(Zind)*1e-6);
            case 2, fprintf('Manual fitting to V(out)(f) normalized near %0.1f MHz...\n',f_data(Zind)*1e-6);
            case 3, fprintf('Manual fitting to Amp(f) normalized near %0.1f MHz...\n',f_data(Zind)*1e-6);
            case 4, fprintf('Manual fitting to phase(f) normalized near %0.1f MHz...\n',f_data(Zind)*1e-6);
            otherwise fprintf('Manual fitting to ratio -V(in)/V(out)...\n');
        end
    else
        fprintf('Manual fitting beam offset V(out) data (normalized to unity)...\n');
    end

    [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj,N,V] = ...
              FDTR_MANUALFIT_vH4(XguessIJ, datparams,...
            sysparams, calparams, matparams, Tparams);

    % Update thermal model in MAIN workspace, so user can access the final
    % phi_model, ratio_model, etc. Only works for LCTE fits.
    XguessIJ(:,1) = Xsol;
    for i = 1:length(XguessIJ(:,1))
        LCTE(XguessIJ(i,2),XguessIJ(i,3)) = XguessIJ(i,1);
    end
    matparams{1} = LCTE;
    
    FDTR_MAIN_MODEL_vH4;
else
    if doughnut, warning('Sorry, autodoughnut is not coded yet.'); break; end
    Xguess = XguessIJ(:,1);

    fprintf('Please wait for automatic fitting...\n');
    tic
    Xsol = fminsearch(@(X) FDTR_FIT_vH4(X, Xij, datparams,...
                           sysparams, calparams, matparams, Tparams),...
                           Xguess,options);

    % fminsearch doesn't output anything but Xsol, so get the rest here.                
    [Z,deltaR_model,ratio_model,LCTE,T_adj]=...
        FDTR_FIT_vH4(Xsol, Xij,datparams,sysparams,calparams,matparams,Tparams);
    toc
end
fprintf('Data fit completed\n');