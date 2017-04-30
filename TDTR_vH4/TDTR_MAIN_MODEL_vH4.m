%% Subfunction of TDTR_MAIN; generates models for time delay, offset, deltaR, and ratio.
% Also vectorizes r_pump based on psc variable. 
%Builds sysparams, calparams, datparams cells.
r_pump = r_probe;
sysparams = {tau_rep f r_pump r_probe}; % initialize sysparams

if ~doughnut
    %vector of time delays (used to generate sensitivity plots)
    if ~exist('N_tdelay','var'), N_tdelay = 20; end
    
    tdelay=logspace(log10(tdelay_min),log10(tdelay_max),N_tdelay)';
    tdelay=linspace(tdelay_min,tdelay_max,N_tdelay)';
    
    if spinning, offset = velocity * tdelay * 1e6; else offset = 0; end % velocity is defined in analyze script, isn't used in FIT scripts.
else
    Nx = 10; % number of offset points to one side of the peak
    tdelay = -20; % in picoseconds; beam offset is generally done at negative delay time.
    %vector of beam offset positions
    offset_span = ceil(2e6 * sqrt(r_probe*r_pump)); % 2x the spot size, rounded up, in microns
    offset = linspace(0, offset_span, Nx)';
end

if psc && length(tdelay) > 1 
    % second condition fixes r_pump at initial value if just 1 time delay.
    % psc: TRUE if pump spot changing over time delay. Relies on
    % "frac", the fractional change, assigned in analyze script.
    r_pump = r_pump*(1 + frac*tdelay/tdelay(end));
    sysparams = {tau_rep f r_pump r_probe}; % update sysparams
end

% Initial and/or normalization time delay, controlled by Zdelay.
[~,Zind] = min(abs(tdelay - Zdelay*1e-12));
calparams = {Zind sigfit intscheme nnodes consider_error};

% Generate initial signal models
[deltaR_model,ratio_model]=TDTR_REFL_vH4(tdelay,matparams,sysparams,A_pump,intscheme,nnodes, offset);
%[deltaR_model_neg,ratio_model_neg]=TDTR_REFL_vH4(-100e-12,matparams,sysparams,A_pump,intscheme,nnodes, 0);

if doughnut % beam offset
    % deltaR and ratio may be a matrix (tdelay, offset). Pull out
    % first (only) time delay element, transpose into column.
    Sig_Amp = deltaR_model(1,:)'; 
    Sig_ratio = ratio_model(1,:)';

    %reconstruct the rest by symmetry
    offset = [-flipud(offset(2:end));offset]; % now offset covers full +/- offset
    deltaR_model = [flipud(Sig_Amp(2:end));Sig_Amp];
    ratio_model = [flipud(Sig_ratio(2:end));Sig_ratio];
    Vout_model = imag(deltaR_model);
    Vin_model = real(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);
    
    
    datparams = {tdelay Vout_model datadir offset};
else % conventional TDTR
    Vin_model = real(deltaR_model);% - real(deltaR_model_neg);
    Vout_model = imag(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);% - sqrt(real(deltaR_model_neg).^2 + imag(deltaR_model_neg).^2);
    
    if ~exist('datadir','var'), datadir = pwd; end
    switch sigfit
        case 1, datparams = {tdelay Vin_model datadir offset};
        case 2, datparams = {tdelay Vout_model datadir offset};
        case 3, datparams = {tdelay Amp_model datadir offset};
        otherwise datparams = {tdelay ratio_model datadir offset}; 
    end
end