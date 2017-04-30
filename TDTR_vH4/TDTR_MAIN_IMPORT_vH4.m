%% Subfunction of TDTR_MAIN, for importing processed TDTR data.
% Input: tdelay, datadir, datain, mid, sigfit, 
%        psc, doughnut, r_pump scalar, Zdelay,
%        calparams, Vout_offset*, detV_raw*, 
%        Voutlinfit*, fitOK*, (*:optional).
%
% Output: signal_data and offset arrays,
%         r_pump array, datparams, calparams.

DM1 = dlmread(datain);
if SJ2
    offset = 0; % no beam offset at SJ2
    warning('V(in) and V(out) offset are not applied for amplitude data set from SJ2');
    if doughnut, error('This code is not prepared for beam offset from SJ2 data'); end
    
    tdelay_raw = DM1(1:50:end,1)*1e-12;
    Amp_raw = DM1(1:50:end,2);
    
    [tdelay,Amp_data] = extract_interior(tdelay_raw,Amp_raw,tdelay_min,tdelay_max);
    [~,Zind] = min(abs(tdelay - Zdelay*1e-12));
    
    Vin_data = Amp_data; % I don't know if SJ2 is R or per-pulse V(in), because sample is spinning.
    
    if psc % psc: is pump spot changing over time delay? Relies on
           % "frac", the fractional change, assigned in analyze script.
        r_pump = r_pump*(1 + frac*tdelay/tdelay(end));
        sysparams = {tau_rep f r_pump r_probe}; % update sysparams
    end
else % not SJ2 data format
    if ~doughnut % conventional TDTR
        offset = 0; % pump and probe beams are overlapped.

        tdelay_raw=DM1(:,2)*1e-12; %imported in picoseconds.  Change to seconds.
        Vin_raw=DM1(:,3); 
        Vout_raw=DM1(:,4);
        ratio_raw=DM1(:,5);

        [~,ratio_data]          =extract_interior(tdelay_raw,ratio_raw, tdelay_min,tdelay_max);
        [~,Vin_data]            =extract_interior(tdelay_raw,Vin_raw,   tdelay_min,tdelay_max);
        [tdelay,Vout_data] =extract_interior(tdelay_raw,Vout_raw,  tdelay_min,tdelay_max);
        
        Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
        
        if psc % psc: is pump spot changing over time delay? Relies on
               % "frac", the fractional change, assigned in analyze script.
            r_pump = r_pump*(1 + frac*tdelay/tdelay(end));
            sysparams = {tau_rep f r_pump r_probe}; % update sysparams
        end

        

        % If the analysis comes with a reasonable linear fit to V_out, 
        % use that to smooth ratio_data. [NOT TESTED]
        if exist('Voutlinfit','var') && exist('fitOK','var')
            if Voutlinfit && fitOK
                Vout_lin = m1*1e12 * tdelay + m2; % [m1 uV/ps]*[1e12 ps/s]*[tdelay seconds] + uV
                r_lin = -Vin_data ./ Vout_lin;
                ratio_data = r_lin;
            end
        end
    else % beam offset TDTR

        offset = DM1(:,1);
        Vout_data = DM1(:,4);
        Vin_data = DM1(:,3);
        ratio_data = -Vin_data./Vout_data;

        offset = offset-mid; % mid is defined in the process or analyze script
    end
    
    % Define Zind: Zdelay = tdelay(Zind), approximately.
    [~,Zind] = min(abs(tdelay - Zdelay*1e-12));

    % shifting V(out) up or down to account for false signal component
    if exist('Vout_offset','var')
        Vout_data = Vout_data - Vout_offset;
        ratio_data = -Vin_data ./ Vout_data;
        Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
    end

    % shifting V(in) up or down to account for false signal component
    if exist('Vin_offset','var')
        Vin_data = Vin_data - Vin_offset;
        ratio_data = -Vin_data ./ Vout_data;
        Amp_data = sqrt(Vin_data.^2 + Vout_data.^2);
    end % (added Feb. 7th, 2015)

    % introducing detector voltage for dR/dT and such.
    if length(DM1(1,:)) > 5 && ~doughnut
        detV_raw = DM1(:,6); % correct for beam offset data, too.
    end 
    if exist('detV_raw','var') && ~doughnut
        [~,detV_data] =extract_interior(tdelay_raw,detV_raw,tdelay_min,tdelay_max);
    end
end

if spinning 
    % now that I have the data inputted, is the sample spinning? 
    % If so, I need to assign an offset vector as a function of tdelay.
    offset = velocity * tdelay * 1e6; % generate offset in MICRONS, for consistency with beam offset input data format.
end
    
%% Compose remaining parameter cells, 
% now that Zind, tdelay, and data are defined.
switch sigfit
    case 1, datparams = {tdelay Vin_data datadir offset};
    case 2, datparams = {tdelay Vout_data datadir offset};
    case 3, datparams = {tdelay Amp_data datadir offset};
    otherwise datparams = {tdelay ratio_data datadir offset}; 
end

calparams = {Zind sigfit intscheme nnodes consider_error};