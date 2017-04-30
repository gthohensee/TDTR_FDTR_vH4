
%% Define the directory and filenames
directory = 'pwd'; % pwd is the current directory
yymmdd = '160429'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('\',yymmdd,'\'); % folder in directory where processed data is kept
savefolder = strcat('\',yymmdd,'_edit\'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

tagline = 'Au'; % a shared prefix for your series of TDTR data files.

filenames = dir(strcat(savedir,tagline,'*'));
nfiles = length(filenames); % number of data files
fstr = '';
for qq = 1:nfiles
    fstr = strcat(fstr, filenames(qq).name,'\n');
end
sprintf(fstr)
%% General modeling conditions
FDTR_init_default; % not BI
P = -1;
P0 = -1;
doughnut = 1;

f_min = 0.02e6;
f_max = 25e6;
N_f = 40;

%% Arrays for dataset-dependent information

thickness = [83.835 2.257 1000 500e3
             89.769 3.182 300 500e3
             87.835 2.257 300 500e3];
stack = {   'Au24K','Ta','SiO2','Si';
            'Au24K','Ta','SiO2','Si';
            'Au24K','Ta','SiO2','Si'};

        
% each column of aniso is TRUE if corresponding layer is anisotropic
aniso = zeros(1,length(thickness(1,:))); 

% A 1D array of data sets requires 1D arrays of measurement parameters.
r_list = 20 * ones(1,nfiles);
jabs_list   = 0 * ones(1,nfiles); % column index of absorption layer
jtrans_list = jabs_list + 1; % column index of transducer layer
measured_pump = 100e-3 * ones(1,nfiles);  % raw powermeter reading, W
measured_probe = 3.7e-3 * ones(1,nfiles); % raw powermeter reading, W

% Xij(m,1:3) = [i j]: assigns LCTE(i,j) as mth fit parameter.
%  Fitting more than two parameters is best done in an exploratory context,
% since TDTR data generally specifies at most two parameters uniquely.
%
%  To fit different parameters for different samples, edit Xij
% and re-run the for loop below for different filename indices.
clear Xij;
Xij(1,1:2) = [1,1]; % [1,3] = thermal conductivity of layer #3
Xij(2,1:2) = [1,2]; % [1,4] = thermal conductivity of layer #4
Xij(3,1:2) = [1,3]; % [3,2] = thickness of layer #2
Xij(4,1:2) = [3,1]; % [3,4] = thickness of layer #4
Xij(5,1:2) = [3,2];
Xij(6,1:2) = [3,3];
Xij(7,1:2) = [3,4];
%Xij(5,1:2) = [3,4];
%Xij(6,1:2) = [3,3];
%% Initialize output; careful not to erase your fit results!
%output = zeros(nfiles,5);
nf = length(Xij(:,1)); % number of fit parameters
Xguess = zeros(1,nf);
XguessIJ = zeros(nf,3);
%% Call prev_output
caltag = '';

switch sigfit
    case 1, caltag = strcat(caltag,'_vinfit');
    case 2, caltag = strcat(caltag,'_voutfit');
    case 3, caltag = strcat(caltag,'_ampfit');
    case 4, caltag = strcat(caltag,'_phasefit');
    otherwise, caltag = strcat(caltag,'_rfit');
end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

prev_output = dlmread(strcat(datadir, yymmdd, '_solutions_',tagline,caltag,'.txt'));
%% Perform fitting.

for ii = 1:3
   
   fname = filenames(ii).name;
   datain = strcat(savedir, fname)
   
   jabs = jabs_list(ii);     % index of absorption layer
   jtrans = jtrans_list(ii); % index of transducer layer
   
   % Laser powers that reach a typical in-air sample.
   tc = 0.8;
   A_pump = measured_pump(ii)*pm*tc;     
   A_probe = measured_probe(ii)*2*pm*tc; % 2x for optical chopper

   %Construct the thermal material properties matrix.
   %refdir is the path to a folder containing all (interpolated, as
   %necessary) reference information or data for building the T-dependent
   %arrays of thermal parameters in T_LCTE. [Not Implemented in FDTR]
   refdir = '';
   
   [LCTE,T_LCTE] = writeLCTE_vH4(stack(ii,:), thickness(ii,:),Xij,P0,T0,refdir);
   % could also have written: LCTE = writeLCTE_vH4(stack,thickness).
   
   fprintf('transducer thickness is %f nm\n', thickness(ii,jtrans_list(ii)));
   
   %Which variable(s) are you fitting for?
   for j = 1:nf
       Xguess(j) = LCTE(Xij(j,1),Xij(j,2));
       XguessIJ(j,:) = [Xguess(j), Xij(j,1), Xij(j,2)];
   end
   
   % Set initial guess for first two fit parameters; useful for P- or T- 
   % series or re-analyzing data. Uncomment the next three lines and set 
   % the appropriate indices for prev_output and LCTE.
   %XguessIJ(1:2,1) = prev_output(ii,2:3);
   %LCTE(1,3:4) = prev_output(ii,2:3);
   %XguessIJ(1,1) = 0.2;
   %LCTE(1,3) = 0.2;
   %LCTE(1,2) = 1;
   
   clear Vin_offset; clear Vout_offset; clear f_data;
   
   sigfit = 4; % 1 == Vin, 2 == Vout, 0 = Ratio, 3 == amp, 4 == phase
   manualfit = 1;
   importdata = 1;
   senseplot = 0;
   doughnut = 0;
   
   switch r_list(ii)
       case 20
           r_pump = 2.0e-6;
           r_probe = r_pump;
       case 5
           r_pump = 8.64e-6;
           r_probe = r_pump; % "effective", correlated spot size, from boffset fits at 8 MHz.
   end
   
   TCR = 2.4e-4;
   
   %phishift = 180; % ARBITRARY phase shift to phi_data only. Not currently implemented.
   %ZF = 8e6;
   
   FDTR_MAIN_vH4;  % Now that all system and material parameters are
                  % initialized, FDTR_MAIN_vH begins the thermal modeling,
                  % fitting, printing of figures, saving files, and
                  % everything else for this particular data file.
                  
   % write output row indicating fit results for this sample.
   % Most of these variables were generated in TDTR_MAIN_vH4.
   output(ii,1:length(Xsol)+3) = [r_pump*1e6, r_probe*1e6, Xsol',Z];
end
%% Write fit results to a text file in the processed data folder.
caltag = '';

switch sigfit
    case 1, caltag = strcat(caltag,'_vinfit');
    case 2, caltag = strcat(caltag,'_voutfit');
    case 3, caltag = strcat(caltag,'_ampfit');
    case 4, caltag = strcat(caltag,'_phasefit');
    otherwise, caltag = strcat(caltag,'_rfit');
end

if manualfit, caltag = strcat(caltag,'_man');
else caltag = strcat(caltag,'_auto'); end

dlmwrite(strcat(savedir, yymmdd,'_solutions_',tagline,caltag,'.txt'),output);

%------------- END CODE --------------