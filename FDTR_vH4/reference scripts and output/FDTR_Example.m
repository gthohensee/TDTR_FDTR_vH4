% Example of how to call FDTR_REFL_VV1
%clear all
%close all

%Set materials properties of each layer
lambda=[180 0.1 35]; %W/m-K, thermal conductivity (layer 1 is top layer, layer N is substrate)
C=[2.43 0.1 3.1]*1e6; %J/m^3-K (Volumetric heat capacity)
t=[100 1 500e3]*1e-9; %m (thickness of each layer)
eta=ones(1,numel(lambda)); %isotropic layers, eta=k-inplane/k-throughplane; (thermal conductivity anisotropy)

%laser properties
r_pump = 6e-6; %1/e2 radius, pump beam
r_probe = 6e-6; %1/e2 radius, probe beam
A_pump = 50e-3; % pump laser power (W)
A_probe = 5e-3; % probe laser power (W)
A_tot = A_pump + A_probe;

%transducer properties
absorbance = 0.1; %fraction of laser energy absorbed by transducer
TCR = 1e-4; %thermoreflectance coefficient (1/K), 1e-4 for aluminum at 785nm.

% calculate steady state temperature rise (i.e. Am I frying the sample?)
kmin=1/(10000*max(r_pump,r_probe));
kmax=1/sqrt(r_pump^2+r_probe^2)*10; %defines effective integration region
[deltaT_model]=rombint_VV1(@(kvectin) FDTR_TEMP_VV1(kvectin,0,lambda,C,t,eta,r_pump,r_probe,absorbance*A_tot),kmin,kmax,1); 
Tss = abs(deltaT_model) % steady state temperature rise, K

%calcuate and plot the frequency response
f = logspace(2,8,200);
[deltaR,ratio,phi]=FDTR_REFL_VV1(f,lambda,C,t,eta,r_pump,r_probe,TCR,A_pump);

semilogx(f,phi)
hold on;
plot(REF(:,1),REF(:,2),'-r')