function [S_LCTE,S_sys] = FDTR_senseplot_offset_vH4(datparams,sysparams, calparams, matparams, Tparams)
%senseplot_offset_vH4 - Calculates sensitivity plots dlogR/dlogX for thermal
%model. The sens_consider booleans can be edited to change which
%sensitivities are calculated. "senseplot_offset" is an expansion of
%"senseplot_vH2" that allows sensitivity calculations for beam
%offset signals, or alternately for the FWHM of the offset signal.
%
% Inputs:
%    datparams - {f_data ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error Mcell_err T0_err P0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS for details.]
%
% Outputs
%    S_LCTE - sensitivities to one-channel overlayer parameters
%    S_sys  - {S_f,S_r_pump,S_r_probe};
%       S_f, S_r_pump, S_r_probe - system parameter sensitivities
%
% Other m-files required: TDTR_REFL_vH4.m
% Subfunctions: none
% MAT-files required: none
%
% See also: errorbars_vH4.m, parametric_senseplot_vH4.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 8-Apr-2014 - made into function, comments updated,
%                                harmonized with TTM version.
%                   14-July-2014 - vH2. No changes since June 11th.
%                   24-Nov-2014  - fixed typo in handling eta sensitivity
%                   27-Jan-2015  - beam offset compatibility. S vectors
%                                  are with respect to time or offset.
%                   1-Feb-2015   - retitled vH3, includes FWHM mode.
%                   17-Feb-2015  - updated header comments.
%                   18-Aug-2015  - amplitude compatible (not in offset)
%                   2-Sep-2015   - FDTR adaptation.
%                   13-Nov-2015  - enable FDTR sensitivities in offset, not
%                                  not just V(out).
%                   11-Jul-2016  - vH4.
%------------- BEGIN CODE --------------
%%
FDTR_INITIALIZE_CELLPARAMS_vH4; % unpacks/cleans cellparams (the five inputs)
fprintf('Generating Sensitivity Plot for all Variables...please be patient\n')

if doughnut
    by_FWHM = -1;
    while by_FWHM == -1
        by_FWHM = input('Enter 1 if you want to study FWHM sensitivities, else 0 to get regular signal sensitivities: ');
        if isempty(by_FWHM), by_FWHM = 0; end
        
        if by_FWHM ~= 1 && by_FWHM ~= 0
            by_FWHM = -1;
            sprintf('Invalid input, try again.\n');
        end
    end
    
    if by_FWHM && sigfit == 4
        warning('I have not figured out FWHM sensitivities for the phase. It does not fit to a Gaussian (more parabolic). Going back to regular sensitivities.');
        by_FWHM = 0;
    end
end

%% Generates a reference model based on anticipated parameters
% Normalize V(in) and V(out) reference model, if fitting by these.
% The commented out lines are not used, unless you want to revise
% the sensitivity calculation to use raw data as its reference.
% I find that messy, not intrinsic or helpful, so sensitivities are
% defined as model perturbations relative to the model.

if ~doughnut % conventional TDTR
    [deltaR_model,ratio_model,phi_model]=FDTR_REFL_vH4(f_data,matparams,sysparams,A_pump,intscheme,nnodes,0);
    
    Vin_model = real(deltaR_model);
    Vout_model = imag(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);
    
    switch sigfit
        case 1 % V(in) fit
            % Construct normalized V(in) model and data,
            % relative to its value at index Zind (from ZF picoseconds, or offset).
            Vin_model_ZF = Vin_model(Zind);
            nVin_model = Vin_model / Vin_model_ZF;

            %Vin_data_ZF = Vin_data(Zind);
            %nVin_data = Vin_data / Vin_data_ZF;
            
        case 2 % V(out) fit
            % Construct normalized V(out) model and data,
            % relative to its mean value near ZF picoseconds.
            Vout_model_ZF = Vout_model(Zind);
            nVout_model = Vout_model / Vout_model_ZF;

            %Vout_data_ZF = mean(Vout_data(Zind:Zind+3));
            %nVout_data = Vout_data / Vout_data_ZF;
            
        case 3
            Amp_model_ZF = Amp_model(Zind);
            nAmp_model = Amp_model / Amp_model_ZF;
            
            %Amp_data_ZF = Amp_data(Zind);
            %nAmp_data = Amp_data / Amp_data_ZF;
            
        case 4
            % phi_model given by REFL.
        otherwise % ratio fit
            % do nothing here
    end
else % beam offset
    if length(f_data) > 1, error('Senseplot does not handle f-dependent offset calculations. See parametric senseplot instead.'); end
    
    % Thermal model is symmetric around 0 offset, so only use
    % the positive offset side of the offset vector to generate model.
    [~,I] = min(abs(offset));
    xvect = offset(I:length(offset));
    [deltaR_model,ratio_model,phi_model]=FDTR_REFL_vH4(f_data,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
    
    % deltaR and ratio may be a matrix (f_data, offset). Pull out
    % first (only) frequency element, transpose into column.
    Sig_Amp = deltaR_model(1,:)'; 
    Sig_ratio = ratio_model(1,:)';
    Sig_phi = phi_model(1,:)';
    
    %reconstruct the rest by symmetry; no data, so no need to interpolate.
    offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
    deltaR_model = [flipud(Sig_Amp(2:end));Sig_Amp];
    
    Vin_model = real(deltaR_model);
    Vout_model = imag(deltaR_model);
    Amp_model = sqrt(Vin_model.^2 + Vout_model.^2); % not to be confused with Sig_Amp, which is deltaR.
    ratio_model = [flipud(Sig_ratio(2:end));Sig_ratio];
    phi_model = [flipud(Sig_phi(2:end));Sig_phi];
    
    
    [FWHM_model, nVin_model, nVout_model, nAmp_model, nphi_model] = do_the_normalize(sigfit, by_FWHM, offset, Vin_model, Vout_model, Amp_model, phi_model, ratio_model);
end

%% initialize matrices for sensitivities
LCTEtemp = LCTE;

%nt = length(f_data); 
nl = length(LCTE(1,:));

S_LCTE = cell(4,nl);

deltaR_temp  = cell(4,nl);
ratio_temp   = cell(4,nl);
phi_temp = cell(4,nl);
%% which sensitivities to consider? (saves time)
LCTE_sens_consider = ones(4,nl);

% Skip absorption layer; instead, couple its perturbations into that of
% the transducer layer. If there is no transducer layer, just a "1 nm"
% absorption layer, set jabs = 0 in your analyze script.
if jabs ~= 0, LCTE_sens_consider(:,jabs) = 0; end 

for j = 1:nl
    % Skip thicknesses and heat capacities of all interfaces in LCTE
    if LCTE(2,j) == 1e5 && LCTE(3,j) == 1e-9 % reliable interface indicators
        LCTE_sens_consider(2,j) = 0;
        LCTE_sens_consider(3,j) = 0;
        LCTE_sens_consider(4,j) = 0; % interfaces are not anisotropic.
    end
    
    % skip eta for all isotropic layers
    if ~aniso(j), LCTE_sens_consider(4,j) = 0; end
end
LCTE_sens_consider

%% -----------------Compute sensitivities for LCTE--------------
for i = 1:4
    for j=1:nl
        % skip conditions; includes the aniso variable.
        if LCTE_sens_consider(i,j) == 0, continue; end

        LCTEtemp(i,j) = LCTE(i,j)*1.01;

        % couplings between model parameters, e.g. absorption layers and
        % anisotropies.
        if jabs ~= 0 && j == jabs, LCTEtemp(i,jtrans) = LCTE(i,jtrans)*1.01; end
        if jabs ~= 0 && j == jtrans, LCTEtemp(i,jabs) = LCTE(i,jabs)*1.01; end
        if i == 1 && aniso(j), LCTEtemp(4,j) = LCTE(4,j)/1.01; end 
        % eta and Lz are coupled; need to hold Lx fixed when varying Lz
        
        % Perturbing eta is assumed to be perturbing Lx, not Lz. %

        % Perform sensitivity calculation with LCTEtemp
        matparams{1} = LCTEtemp;
        
        % Thermal model is symmetric around 0 offset, so only use
        % the positive offset side of the offset vector to generate model.
        [~,I] = min(abs(offset));
        xvect = offset(I:length(offset)); % this will be 0 if conventional TDTR
        [deltaR_temp{i,j},ratio_temp{i,j},phi_temp{i,j}] = FDTR_REFL_vH4(f_data,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
        
        
        if ~doughnut
            Vin_temp = real(deltaR_temp{i,j});
            Vout_temp = imag(deltaR_temp{i,j});
            Amp_temp = sqrt(Vin_temp.^2 + Vout_temp.^2);
            
            switch sigfit
                case 1
                    norm = Vin_temp(Zind);
                    nVin_temp = Vin_temp / norm;
                    Num=log(nVin_temp)-log(nVin_model);
                case 2
                    norm = Vout_temp(Zind);
                    nVout_temp = Vout_temp / norm;
                    Num=log(nVout_temp)-log(nVout_model);
                case 3
                    norm = Amp_temp(Zind);
                    nAmp_temp = Amp_temp / norm;
                    Num=log(nAmp_temp)-log(nAmp_model);
                case 4
                    Num=log(phi_temp{i,j})-log(phi_model);
                otherwise
                    Num=log(ratio_temp{i,j})-log(ratio_model);
            end
        else % beam offset
            % deltaR and ratio may be a matrix (f_data, offset). Pull out
            % first (only) time delay element, transpose into column.
            temp_deltaR = deltaR_temp{i,j};
            temp_ratio = ratio_temp{i,j};
            temp_phi = phi_temp{i,j};
            
            Sig_Amp = temp_deltaR(1,:)'; 
            Sig_ratio = temp_ratio(1,:)';
            Sig_phi = temp_phi(1,:)';
            
            %reconstruct the rest by symmetry; no data, so no need to interpolate.
            offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
            temp_deltaR = [flipud(Sig_Amp(2:end));Sig_Amp];
            
            Vin_temp = real(temp_deltaR);
            Vout_temp = imag(temp_deltaR);
            Amp_temp = sqrt(Vin_temp.^2 + Vout_temp.^2);
            temp_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];
            temp_phi = [flipud(Sig_phi(2:end));Sig_phi];
            
            [FWHM_temp, nVin_temp, nVout_temp, nAmp_temp, temp_nphi] = ...
                do_the_normalize(sigfit, by_FWHM, offset, Vin_temp, Vout_temp, Amp_temp, temp_phi, temp_ratio);
            
            % Alternative mode for depicting beam offset sensitivities:
            % instead of sensitivity to V(out) as a function of offset,
            % condense this information into the scalar change in FWHM
            % of V(out) at specified time delay.
            if by_FWHM && length(offset) > 4
                Num = log(FWHM_temp) - log(FWHM_model); % this is scalar
            else
                if by_FWHM && length(offset) <= 4
                    warning('Cannot fit FWHM from so few offset points. Defaulting back to V(out)(x).')
                    by_FWHM = 0;
                    sigfit = 2;
                end
                
                switch sigfit
                    case 1, Num = log(nVin_temp) - log(nVin_model); % this is array of length(offset).
                    case 2, Num = log(nVout_temp) - log(nVout_model); % this is array of length(offset).
                    case 3, Num = log(nAmp_temp) - log(nAmp_model); % this is array of length(offset).
                    case 4
                        Num = log(temp_nphi) - log(nphi_model); % this is array of length(offset).
                        %Num = log(temp_phi) - log(phi_model); % this is array of length(offset).
                    otherwise, Num = log(temp_ratio) - log(ratio_model); % this is array of length(offset).
                end
            end
        end % end doughnut
        
        Denom=log(LCTEtemp(i,j))-log(LCTE(i,j));
        S_LCTE{i,j}=Num/Denom;

        LCTEtemp = LCTE;
        matparams{1} = LCTE;
        fprintf('Calculated S_LCTE{%i,%i}...\n',i,j);
    end
end
clear LCTEtemp;

%% -----------------Compute sensitivities for system parameters-----------
sys_consider = [1 2];
S_r_pump = []; S_r_probe = [];
for i = sys_consider
    sysbump = 1.01;
    switch i
        case 1, sysparams{1} = r_pump*sysbump;
        case 2, sysparams{2} = r_probe*sysbump;
    end
    
    % Thermal model is symmetric around 0 offset, so only use
    % the positive offset side of the offset vector to generate model.
    [~,I] = min(abs(offset));
    xvect = offset(I:length(offset)); % this will be 0 if conventional TDTR
    [temp_deltaR,temp_ratio,temp_phi]=FDTR_REFL_vH4(f_data,matparams,sysparams,A_pump,intscheme,nnodes,xvect);
    if ~doughnut
        Vin_temp = real(temp_deltaR); 
        Vout_temp = imag(temp_deltaR);
        Amp_temp = sqrt(Vin_temp.^2+Vout_temp.^2);
        switch sigfit
            case 1
                norm = Vin_temp(Zind);
                nVin_temp = Vin_temp / norm;
                Num=log(nVin_temp)-log(nVin_model);
            case 2
                norm = Vout_temp(Zind);
                nVout_temp = Vout_temp / norm;
                Num=log(nVout_temp)-log(nVout_model);
            case 3
                norm = Amp_temp(Zind);
                nAmp_temp = Amp_temp / norm;
                Num=log(nAmp_temp)-log(nAmp_model);
            case 4
                Num=log(temp_phi)-log(phi_model);
            otherwise
                Num=log(temp_ratio)-log(ratio_model);
        end
    else % beam offset, assumes f_data has just one element.
        % deltaR and ratio may be a matrix (f_data, offset). Pull out
        % first (only) time delay element, transpose into column.
        Sig_Amp = temp_deltaR(1,:)'; 
        Sig_ratio = temp_ratio(1,:)';
        Sig_phi = temp_phi(1,:)';

        %reconstruct the rest by symmetry; no data, so no need to interpolate.
        offset = [-flipud(xvect(2:end));xvect]; % makes sure offset and model line up.
        temp_deltaR = [flipud(Sig_Amp(2:end));Sig_Amp];
    
        Vin_temp = real(temp_deltaR);
        Vout_temp = imag(temp_deltaR);
        Amp_temp = sqrt(Vin_temp.^2 + Vout_temp.^2);
        temp_phi = [flipud(Sig_phi(2:end));Sig_phi];
        temp_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];
        
        [FWHM_temp, nVin_temp, nVout_temp, nAmp_temp, temp_nphi] = ...
                do_the_normalize(sigfit, by_FWHM, offset, Vin_temp, Vout_temp, Amp_temp, temp_phi, temp_ratio);
            
        % Alternative mode for depicting beam offset sensitivities:
        % instead of sensitivity to V(out) as a function of offset,
        % condense this information into the scalar change in FWHM
        % of V(out) at specified time delay.
        if by_FWHM && length(offset) > 4
            Num = log(FWHM_temp) - log(FWHM_model); % this is scalar
        else
            if by_FWHM && length(offset) <= 4
                warning('Cannot fit FWHM from so few offset points..')
                by_FWHM = 0;
            end

            switch sigfit
                case 1, Num = log(nVin_temp) - log(nVin_model); % this is array of length(offset).
                case 2, Num = log(nVout_temp) - log(nVout_model); % this is array of length(offset).
                case 3, Num = log(nAmp_temp) - log(nAmp_model); % this is array of length(offset).
                case 4
                    figure(1)
                    clf
                    plot(offset, temp_nphi, 'r-', offset, nphi_model,'k-')
                    Num = log(temp_nphi) - log(nphi_model); % this is array of length(offset).
                    %Num = log(temp_phi) - log(phi_model); % this is array of length(offset).
                otherwise, Num = log(temp_ratio) - log(ratio_model); % this is array of length(offset).
            end
        end
        
    end
    Denom=log(sysbump);
    
    switch i
        case 1, S_r_pump=Num/Denom;
        case 2, S_r_probe=Num/Denom;
    end
    sysparams = {r_pump, r_probe}; % return to reference value
    fprintf('Calculated S_sys #%i...\n',i);
end

%% Plot sensitivities
figure(203)
clf
is_fscan = ~doughnut || by_FWHM;

if is_fscan, axes('XScale','log'); end
set(gca,'Box','on');
hold on;
ColorOrder = get(gcf,'DefaultAxesColorOrder');

%% label LCTE sensitivities
LCTElegend = [];
LCTElab = cell(size(LCTE));
LCTEmarker = {'o','*','x','+'};

for i = 1:4
    for j = 1:nl
        switch i
            case 1, LCTElab{i,j} = sprintf('L%i',j);
            case 2, LCTElab{i,j} = sprintf('C%i',j);
            case 3, LCTElab{i,j} = sprintf('t%i',j);
            case 4, LCTElab{i,j} = sprintf('e%i',j);
        end
        
        if isempty(S_LCTE{i,j}), continue;
        else
            if ~doughnut || by_FWHM,
                if sigfit == 4 % negative signals behave funny in sensitivity calculation; flip sign. 
                               % Normalized signals are all positive. What
                               % about ratio in FDTR?
                    semilogx(f_data/1e6,-1*S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
                else
                    semilogx(f_data/1e6,-1*S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
                end
            else
                plot(offset,S_LCTE{i,j},LCTEmarker{i},'Color',ColorOrder(j,:));
            end
            LCTElegend = [LCTElegend;LCTElab{i,j}];
        end
    end
end

%% plot and label system sensitivities
syslegend = [];
for i = sys_consider
    switch i
        case 1
            if is_fscan,
                semilogx(f_data/1e6,S_r_pump,'-','Color',ColorOrder(i,:));
            else
                plot(offset,S_r_pump,'-','Color',ColorOrder(i,:));
            end
            syslegend = [syslegend;'Rp'];
        case 2
            if is_fscan,
                semilogx(f_data/1e6,S_r_probe,'-','Color',ColorOrder(i,:));
            else
                plot(offset,S_r_probe,'-','Color',ColorOrder(i,:));
            end
            syslegend = [syslegend;'Rb'];
    end
end

%% Other plot details
figure(203);
set(gca,'FontSize',16);
set(gca, 'TickLength' , [.02 .02]);

legend([LCTElegend;syslegend])
if is_fscan,
    xlabel('frequency (MHz)','FontSize',16);
else
    xlabel('offset (microns)','FontSize',16);
end

if is_fscan
    if ~doughnut
        switch sigfit
            case 1, ylabel(sprintf('Sensitivity:  dlog[nV(in) @ %0.0f MHz]/dlogX',f_data(Zind)*1e-6), 'FontSize',16);
            case 2, ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f MHz]/dlogX',f_data(Zind)*1e-6), 'FontSize',16);
            case 3, ylabel(sprintf('Sensitivity:  dlog[nAmp @ %0.0f MHz]/dlogX',f_data(Zind)*1e-6), 'FontSize',16);
            case 4, ylabel(sprintf('Sensitivity: -dlog[phase (degrees)]/dlogX'), 'FontSize',16);
            otherwise, ylabel('Sensitivity:  dlogR/dlogX', 'FontSize',16);
        end
    else % by_FWHM = TRUE
        switch sigfit
            case 1, ylabel(sprintf('Sensitivity:  dlog[V(in) FWHM]/dlogX'), 'FontSize',16);
            case 2, ylabel(sprintf('Sensitivity:  dlog[V(out) FWHM]/dlogX'), 'FontSize',16);
            case 3, ylabel(sprintf('Sensitivity:  dlog[amplitude FWHM]/dlogX'), 'FontSize',16);
            case 4, ylabel(sprintf('Sensitivity: -dlog[phase FWHM]/dlogX'), 'FontSize',16);
            otherwise, ylabel(sprintf('Sensitivity:  dlog[ratio FWHM]/dlogX'), 'FontSize',16);
        end
        
    end
    set(gca, 'XTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50]);
    set(gca, 'XTickLabel', [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]);
    set(gca, 'XMinorTick','off')
    axis([min(f_data)/1e6 2*max(f_data)/1e6 -1 1]);
else % beam offset
    switch sigfit
        case 1, ylabel(sprintf('Sensitivity:  dlog[nV(in) @ %0.0f MHz]/dlogX',f_data(1)*1e-6), 'FontSize',16);
        case 2, ylabel(sprintf('Sensitivity:  dlog[nV(out) @ %0.0f MHz]/dlogX',f_data(1)*1e-6), 'FontSize',16);
        case 3, ylabel(sprintf('Sensitivity:  dlog[nAmp @ %0.0f MHz]/dlogX',f_data(1)*1e-6), 'FontSize',16);
        case 4, ylabel(sprintf('Sensitivity: -dlog[phase (degrees) %0.0f MHz]/dlogX',f_data(1)*1e-6), 'FontSize',16);
        otherwise, ylabel(sprintf('Sensitivity: dlogR/dlogX @ %0.0f MHz]/dlogX',f_data(1)*1e-6), 'FontSize',16);
    end
    axis([min(xvect) max(xvect) -1 1]);
end

%% export sensitivities
%S_LCTE = S_LCTE;

S_sys = {S_r_pump,S_r_probe};
end


%% subfunction for normalization of signals, used frequently.
function [FWHM, nVin, nVout, nAmp, nphi] = do_the_normalize(sigfit, by_FWHM, offset, Vin, Vout, Amp, phi, ratio)

FWHM = 0;
nVin = 0;
nVout = 0;
nAmp = 0;
nphi = 0;

switch sigfit
    case 1 % Vin
        sig = Vin;
    case 2 % Vout
        sig = Vout;
    case 3 % amplitude
        sig = Amp;
    case 4 % phase
        sig = phi;
    otherwise % ratio (0)
        sig = ratio;
end

[~,I] = min(abs(offset)); % get index of offset = 0 peak
nsig = sign(sig(I)) * sig / sig(I); % normalize peak to +1

if by_FWHM % since I defined nsig_model, this can be used for any sigfit
    [~,J] = min(abs(nsig - 0.5)); % gets half-rise index.
    Xguess = 0.4247 * offset(J); % if offset was continuous
                                 % instead of a few points,
                                 % this would be the FWHM.

    if ~exist('options','var'); options = optimset('TolFun',1e-2,'TolX',1e-2); end % don't need high precision in fit
    Xsol = fminsearch(@(X) SpotSize_V4(X,nsig,offset),Xguess,options);
    FWHM = 2 * sqrt(2*log(2)) * Xsol;
end

switch sigfit % back to specific variables
    case 1 % Vin
        nVin = nsig;
    case 2 % Vout
        nVout = nsig;
    case 3 % amplitude
        nAmp = nsig;
    case 4
        nphi = nsig;
    otherwise % ratio isn't supposed to be normalized.
        % do nothing
end

end
