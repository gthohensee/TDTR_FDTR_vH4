%% Type and conditions of modeling
P0 = -1;         % GPa initial. P = -1 assumes atmospheric pressure.

BI = 0;         % boolean: TRUE if using bidirectional heat flow model
n_toplayer = 3; % # of model layers above the point where heat is deposited.

%% Common sysparams - FDTR system settings
TCR=-2.5e-4; % default coefficient of thermal reflectance (Au)
pm = 1.00; % optical power is ?.??x the reading of the power meter

%% calparams - calculation parameters
ZF = .1e6; % MHz starting time. In automatic fitting, goodness-of-fit
             % is calculated between ZF and f_max. In V(in)
             % fitting, ZF indicates the frequency at which to
             % normalize the V(in) signal.

sigfit = 0; % TRUE if fitting by the V(in) signal instead of ratio.

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
nnodes = 35;  % number of nodes for Legendre-Gauss or Simpson integration;
              % affects numerical accuracy. Don't go below 35 nodes
              % unless you know what you're doing.

%% Tparams - parameters for temperature-dependent samples
T0 = -1; % nominal (K). Set T0 = -1 to ignore laser heating and assume room temperature.

% optional parameters:
absC = [295 0.13]; % Transducer absorption coefficient (0.13 for room temperature Al)
perpulse = 0; % TRUE if accounting for per-pulse heating as well as steady-state

%% other operational parameters for MAIN
% Time delay boundaries for thermal modeling AND data presentation.
f_min= 0.02e6; 
f_max= 25e6;
N_f = 50;

options = optimset('TolFun',1e-2,'TolX',1e-2); % tolerances for autofitting
                                               % and for error bars.


senseplot = 0; %Generate Sensitivity Plot? This option is available
               %dynamically in MANUALFIT.
              
ebar = 0; %Calculate Errorbars? 0 for no, 1 for yes (takes longer)

importdata = 1; % TRUE if fitting data. FALSE if just running sensitivity
                % plots or error bar calculations.
                
manualfit = 1;  % TRUE if fitting manually, FALSE if auto-fitting.

get_dRdT = 0; % TRUE if extracting dR/dT from the TDTR data.
              % Requires accurate measurement of incident pump
              % and probe powers, as well as photodiode voltage vs. delay
              % time. Keep this off if you don't fully understand it.