%% Subfunction of FDTR_MAIN, for importing processed TDTR data.
% Input: datadir, datain, mid, sigfit, doughnut, ZF, calparams, 
%        Vout_offset*, Vin_offset*, detV_raw*, 
%        (*:optional).
%
% Output: signal_data and offset arrays,
%         datparams, calparams.

DM1 = dlmread(datain);

if ~doughnut % conventional FDTR
    offset = 0; % pump and probe beams are overlapped.

    f_raw=DM1(:,1); %imported in Hz.
    Vin_raw=DM1(:,2);
    Vout_raw=DM1(:,3);
    Amp_raw=DM1(:,4); % amplitude in volts
    phi_raw=DM1(:,5); % phase in degrees

    [f_data,Vin_data]       =extract_interior(f_raw,Vin_raw,   f_min,f_max);
    [~,Vout_data]           =extract_interior(f_raw,Vout_raw,  f_min,f_max);
    [~,Amp_data]            =extract_interior(f_raw,Amp_raw,   f_min,f_max);
    [~,phi_data]          =extract_interior(f_raw, phi_raw, f_min,f_max);
    
    % append nonsense data points so that sensitivities will go to 25 MHz.
   %  f_data(end+1:end+1+10) = [9 10 11 12 14 16 18 20 22 24 25]'*1e6;
   %  Vin_data(end+1:end+1+10) = 1;
   %  Vout_data(end+1:end+1+10) = 0;
   %  Amp_data(end+1:end+1+10) = 1;
   %  phi_data(end+1:end+1+10) = 0;
    
    jj = sqrt(-1);
    deltaR_data = Vin_data + jj*Vout_data;
    ratio_data = -Vin_data ./ Vout_data;
    %phi_data = -angle(deltaR_data)*(180/pi); %phase delay, degrees
    if exist('Amp_offset','var')
        Amp_data = Amp_data - Amp_offset;
    end
    
    % Define Zind: ZF = tdelay(Zind), approximately.
    [~,Zind] = min(abs(f_data - ZF));
else % beam offset FDTR
    offset = DM1(:,1)*.95;
    Vin_data = DM1(:,3); % (:,2) is photodiode voltage
    Vout_data = DM1(:,4);
    Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
    
    
    ratio_data = DM1(:,5);
    
    deltaR_data = Vin_data + sqrt(-1)*Vout_data;
    phi_data = angle(deltaR_data*exp(sqrt(-1)*phishift))*(180/pi); %phase delay, degrees
    
    offset = offset-mid; % mid is defined in the process or analyze script
    Zind = 0;
end

% shifting V(out) up or down to account for false signal component
if exist('Vout_offset','var')
    Vout_data = Vout_data - Vout_offset;
    ratio_data = -Vin_data ./ Vout_data;
    Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
    
    deltaR_data = Vin_data + ii*Vout_data;
    phi_data = -angle(deltaR_model)*(180/pi); %phase delay, degrees
end

% shifting V(in) up or down to account for false signal component
if exist('Vin_offset','var')
    Vin_data = Vin_data - Vin_offset;
    ratio_data = -Vin_data ./ Vout_data;
    Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
    
    deltaR_data = Vin_data + ii*Vout_data;
    phi_data = -angle(deltaR_model)*(180/pi); %phase delay, degrees
end % (added Feb. 7th, 2015)

% introducing detector voltage for dR/dT and such.
clear detV_raw;
if length(DM1(1,:)) > 5 && ~doughnut
    detV_raw = DM1(:,6); % correct for beam offset data, too.
end 
if exist('detV_raw','var') && ~doughnut
    [~,detV_data] =extract_interior(f_raw,detV_raw,f_min,f_max);
end

%% Compose remaining parameter cells, 
% now that Zind, tdelay, and data are defined.
switch sigfit
    case 1, datparams = {f_data Vin_data datadir offset};
    case 2, datparams = {f_data Vout_data datadir offset};
    case 3, datparams = {f_data Amp_data datadir offset};
    case 4, datparams = {f_data phi_data datadir offset};
    otherwise, datparams = {f_data ratio_data datadir offset}; 
end

calparams = {Zind sigfit intscheme nnodes consider_error};