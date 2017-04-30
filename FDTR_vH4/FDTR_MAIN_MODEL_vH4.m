%% Subfunction of FDTR_MAIN; generates models for f, offset, deltaR, and ratio.
%Builds sysparams, calparams, datparams cells.

if ~exist('f_data','var')
    %vector of frequencies (used to generate sensitivity plots)
    f_data = logspace(log10(f_min),log10(f_max),N_f)';
end
if ~exist('datadir','var'), datadir = pwd; end

sysparams = {r_pump r_probe}; % initialize sysparams

if ~importdata
    clear offset;
end

if ~doughnut
    offset = 0;
else
    if ~exist('offset','var')
        flip = 1;
        Nx = 100; % number of offset points to one side of the peak
        %vector of beam offset positions
        offset_span = 0.4;%ceil(3e6*sqrt(r_pump(1)*r_probe)); % Normal FDTR: 3x the geometric mean spot size, rounded up, in microns
        offset = linspace(0, offset_span, Nx)';
    else
        flip = 0;
    end
end

% Initial and/or normalization frequency, controlled by ZF.
[~,Zind] = min(abs(f_data - ZF));
calparams = {Zind sigfit intscheme nnodes consider_error};

% Generate initial signal models
[deltaR_model,ratio_model,phi_model]=FDTR_REFL_vH4(f_data,matparams,sysparams,A_pump,intscheme,nnodes, offset);

if doughnut % beam offset
    % deltaR and ratio may be a matrix (tdelay, offset). Pull out
    % first (only) time delay element, transpose into column.
    Sig_Amp = deltaR_model(1,:)'; 
    Sig_ratio = ratio_model(1,:)';

    if flip
        %reconstruct the rest by symmetry
        offset = [-flipud(offset(2:end));offset]; % now offset covers full +/- offset
        deltaR_model = [flipud(Sig_Amp(2:end));Sig_Amp];
        ratio_model = [flipud(Sig_ratio(2:end));Sig_ratio];
    end
    Vout_model = imag(deltaR_model);
    Vin_model = real(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);
    phi_model = -angle(deltaR_model)*(180/pi); %phase delay, degrees, already flipped for symmetry by deltaR_model.
    
    datparams = {f_data Vout_model datadir offset};
else % conventional FDTR
    Vin_model = real(deltaR_model);
    Vout_model = imag(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);
end
% 
% % Plot amplitude model
% figure(15)
% semilogy(offset, Amp_model/max(Amp_model),'-');
% axis([min(offset) max(offset) -0.2 1.4]);
% xlabel('x (microns)','FontSize',20)
% ylabel('Offset normalized amplitude','FontSize',20)
% hold on;

switch sigfit
    case 1, datparams = {f_data Vin_model datadir offset};
    case 2, datparams = {f_data Vout_model datadir offset};
    case 3, datparams = {f_data Amp_model datadir offset};
    case 4, datparams = {f_data phi_model datadir offset};
    otherwise, datparams = {f_data ratio_model datadir offset}; 
end