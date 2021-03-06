%FDTR_sense_setup - Template for calling FDTR_parametric_senseplot_vH3.
%
% Inputs:
%    none
%
% Outputs
%    See outputs of parametric_senseplot_vH3.
%
% Other m-files required: parametric_senseplot_vH3.m
% Subfunctions: none
% MAT-files required: none
%
% See also: sense_setup_retro.m, parametric_senseplot_vH3.m

% Author: Gregory Hohensee
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 17-Feb-2015 - header comments added.
%                   2-Sep-2015  - FDTR adaptation. No tau_rep, etc.
%------------- BEGIN CODE --------------
%% Define the directory
directory = '\Users\hohensee_g\Documents\MATLAB\TDTR_vH3'; % pwd is the current directory
yymmdd = '150902'; % e.g., 140227 = February 27th, 2014.
savefolder = strcat('\',yymmdd,'_parametric\'); % folder in directory where sensitivities are saved
savedir = strcat(directory, savefolder);

%% Beam intensity and spot size conditions
objective = 20; % strength of objective lens

switch objective % set pump 1/e^2 intensity focused radius, in meters.
    % 12/4/2014: measured 11.5 - 12.5 micron ud - lr spot size for 5x.
    % "tc": transmission coefficient for this objective
    case 2, r_pump = 2*11.7e-6; tc = 0.87; 
    case 5, r_pump = 10e-6;     tc = 0.90; 
    case 10, r_pump = 6.1e-6;   tc = 0.80;
    case 20, r_pump = 2.65e-6;  tc = 0.70;
    case 50, r_pump = 1.4e-6;  tc = 0.70; % tc = ?? for 50x
    %case 100, r_pump = ?e-6;   tc = 0.??;
end 
r_probe = r_pump;         % probe 1/e^2 radius, meters.

% Laser powers that reach a typical in-air sample.
pm = 1; % deviation of powermeter reading from actual average laser power pre-objective.
A_pump = 10 * 1e-3*pm*tc;
A_probe = 5 * 1e-3*pm*tc;

%% General thermal modeling conditions
ZF = 1e6;     % Normalization frequency (Hz) for V(in) and V(out).
              % Irrelevant if calculating by ratio or beam offset.
              
BI = 0;         % boolean: TRUE if using bidirectional heat flow model.
n_toplayer = 3; % for BI: # of model layers above the point where heat is deposited.

intscheme = 0; % integration scheme: 0 = Legendre-Gauss, 1 = Rombint,
               % 2 = Simpson integration.
               % Beam offset only works with Legendre-Gauss for now.
nnodes = 35;  % number of nodes for Legendre-Gauss or Simpson integration;
              % affects numerical accuracy. Don't go below 35 nodes
              % unless you know what you're doing.

T0 = -1; % nominal (K). Set T0 = -1 to ignore laser heating and assume room temperature.
perpulse = 0; % TRUE if accounting for per-pulse heating as well as steady-state
% optional parameters:
absC = [295 0.13]; % Transducer absorption coefficient (0.13 for room temperature Al)
                   % (used in place of TCR if self-correcting for
                   % temperature-sensitive thermal parameters).

consider_error = 0; % Boolean; 0 if not performing errorbar calculations.
%% other operational parameters for MAIN
TCR=1e-4; % default coefficient of thermal reflectance

% Frequency boundaries for thermal modeling AND data presentation.
f_min= 0.05e6; % warning: the smaller this is, the longer the
                    % computation time of the thermal model.
f_max= 25e6; % approx. max extent of our delay stage.
N_f = 50; % number of logspace time delay points between min and max.
f_data=logspace(log10(f_min),log10(f_max),20); % pump laser modulation frequency, Hz

options = optimset('TolFun',1e-2,'TolX',1e-2); % tolerances for autofitting
                                               % and for error bars.

%% Construct the thermal material properties matrix.
thickness = [100 1 500e3]; % thickness in nm of each model layer
stack = {'SiGe','/','SiO2'}; % identifier for each model layer
Xijc = 0; % LCTE(i,j) indices of the "c" fit parameters; 0 if none.
P0 = -1; % Pressure (GPa); -1 if not high pressure.
T0 = -1; % Temperature (K); -1 if room temperature.
refdir = ''; % directory for T-dependent thermal parameter data files
rho_Al = 4.8e-8; % Al transducer resistivity; -1 for default value.

[LCTE,T_LCTE] = writeLCTE_vH3(stack, thickness,Xijc,P0,T0,refdir,rho_Al);
% could also have written: LCTE = writeLCTE_vH3(stack,thickness).

jabs = 0;           % index of absorption layer (for T-dependence and coupled sensitivities)
jtrans = jabs + 1;  % index of transducer layer (for T-dependence and coupled sensitivities)

% each column of aniso is TRUE if corresponding layer is anisotropic.
% aniso is used in sensitivity calculation to decide how to decouple eta and lambda.
aniso = zeros(1,length(LCTE(1,:)));

%% Choose type of signal: V(out) offset, V(in)(t), V(out)(t), or ratio(t)?
doughnut = 1;   % boolean: TRUE if calculating beam offset sensitivities.
sigfit = 2; % only relevant if not doughnut: 1 == Vin, 2 == Vout, 0 = Ratio
spinning = 0;

%% Assemble cell params.
matparams = {LCTE aniso BI n_toplayer TCR doughnut spinning};
sysparams = {r_pump r_probe};
Tparams = {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans};
% calparams assignment needs to wait until Zind is defined.
% datparams assignment needs to wait until tdelay and data are defined.
% If psc == TRUE (varying spot size), sysparams will be updated later,
%   after tdelay is defined.

%% Tweak LCTE
%LCTE(1,3) = 0.2; % proto 200 MW
%LCTE(2,4) = 2e6;
matparams{1} = LCTE;

%% produce model tdelay and/or offset arrays, plus r_pump array (for psc),
% and pack up the rest of the cell params.
f_data = .1e6;
FDTR_MAIN_MODEL_vH3; 
%% Define independent variable
xvar = logspace(log10(1),log10(100),15); % range of xvar parameter values in LCTE (SI) units
Xij = 0;
Xij = [1,1]; % index of LCTE matrix of xvar parameter.
SYS_i = 0; % Scalar. [1,2,3,4] indicates xvar is [f, r_pump, r_probe, tau_rep], respectively.
           % Zero values in Xij or SYS_i indicate that xvar is not in LCTE
           % or a system parameter, respectively.

%% Which sensitivities to calculate?
LCTE_sens_consider = zeros(4,length(LCTE(1,:)));
LCTE_sens_consider(1,1) = 1;
LCTE_sens_consider(1,3) = 1;
LCTE_sens_consider(2,1) = 1;
sys_consider = [0 0 0 0]; %[r_pump, r_probe, w0, phase]

%% Execute sensitivity calculations and plotting.
tic
[S_LCTE,S_sys,xvar,SS_LCTE,SS_sys] = ...
    FDTR_parametric_senseplot_vH3(Xij,SYS_i, xvar, ...
                   datparams,sysparams, calparams, matparams, Tparams, ...
                   LCTE_sens_consider, sys_consider);
toc
