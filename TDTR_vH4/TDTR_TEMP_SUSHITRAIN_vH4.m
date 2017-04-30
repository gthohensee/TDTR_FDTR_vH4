function Integrand=TDTR_TEMP_SUSHITRAIN_vH4(kvectin,freq,LCTE,r_pump,r_probe,A_pump,xoffset)
%TDTR_TEMP_SUSHITRAIN_vH4 - For beam offset measurements, computes frequency 
% domain average temperature response to periodic gaussian pump beam, 
% probed by another gaussian beam.
%
% This program is vectorized/optimized to handle all frequencies 
% simultaneously (f is a ROW vector)
%
% Syntax:  [Integrand,G]=TDTR_TEMP_vH4(kvectin,freq,LCTE,r_pump,r_probe,A_pump)
%
% Inputs:
%    kvectin - vector of wavenumber (m^-1)
%    freq    - excitation frequency (Hz), ROW vector
%    LCTE    - vertcat(lambda,C,t,eta)
%    lambda  - vector of thermal conductivities, 
%              lambda(1)=top surface,(W/m-K)
%    C       - vector of volumetric specific heat (J/m3-K)
%    t       - thicknesses of each layer (layer N will NOT be used, semiinfinite)
%    eta     - anisotropy of each layer.
%    r_pump  - Pump spot size (m)
%    r_probe - Probe spot size (m)
%    A_pump  - Pump power (W), used to ESTIMATE amplitude (not used for fitting)
%    xoffset - lateral relative offset of pump and probe beams (vector)
%
% Outputs:
%    Integrand - G.*Kernal; see David Cahill's 2004 paper.
%    G         - The layer G(k), where k is spatial wavenumber
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_REFL_vH4.m

% Author: Joseph P. Feser
% --
% email: --
% Website: TDTR_V4 package: http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html
% September 2012; 
% Revision history: 13-November-2014: Merged with vH2 package. Renamed
%                                     EXTRA to xoffset.
%                   17-Feb-2014: vH3. No change.
%                   11-Jul-2016: SUSHI TRAIN is my (failed) attempt at emulating
%                                a TDTR system where the sample was moving at
%                                some velocity (spinning disk). vH4.
%------------- BEGIN CODE --------------
%% Unpacked from my other data structures so I don't need to edit the code.
lambda = LCTE(1,:);
C = LCTE(2,:);
t = LCTE(3,:);
eta = LCTE(4,:);

%% Joe's code from here onward.
Nfreq=length(freq);
kvect=kvectin(:)*ones(1,Nfreq);
Nlayers=length(lambda); %# of layers
Nint=length(kvectin); %# of different frequencies to calculate for

%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

%kmax=1/sqrt(r_pump^2+r_probe^2)*1.5; %cutoff wavevector for integration
ii=sqrt(-1);
alpha=lambda./C;
omega=2*pi*freq;
q2=ones(Nint,1)*(ii*omega./alpha(Nlayers));
kvect2=kvect.^2;

un=sqrt(4*pi^2*eta(Nlayers)*kvect2+q2);
gamman=lambda(Nlayers)*un;
Bplus=zeros(Nint,Nfreq);
Bminus=ones(Nint,Nfreq);
kterm2=4*pi^2*kvect2;
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=ones(Nint,1)*(ii*omega./alpha(n-1));
        unminus=sqrt(eta(n-1)*kterm2+q2);
        gammanminus=lambda(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*Bplus+BB.*Bminus;
        temp2=BB.*Bplus+AA.*Bminus;
        expterm=exp(unminus*t(n-1));
        Bplus=(0.5./(gammanminus.*expterm)).*temp1;
        Bminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(t(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        un=unminus;
        gamman=gammanminus;
    end
end

G=(Bplus+Bminus)./(Bminus-Bplus)./gamman; %The layer G(k)

%% GTH edit: varying pump size over delay time, August 5th 2014
% extend [k,f] matrices G and Kernal into 3D by projecting
% across r_pump, which is a function of time delay.
% Take guidance from this command: reshape(test(:) * [1 2], 2, 2, [])
% The original calculation, for scalar r_pump, was...
% %Kernal=2*pi*P.*kvect; %The rest of the integrand
% %Integrand=G.*Kernal;

% SUSHITRAIN is only called when spinning and not doughnut. Therefore
% xoffset has length tdelay, and the Integrand should have a dimension
% along the tdelay axis. So Nt = length(time delay).

[Nk,Nf] = size(kvect); Nt = length(r_pump); % define dimensions
arg1 = -pi^2*(r_pump.^2)/2; % column vector C, size(Nt,1); no r_probe because r_probe is visualized as the offset beam. Symmetric under pump/probe swap.
arg2 = kvect.^2; % matrix B, size(Nk,Nf)

expterm = exp(reshape(arg2(:) * arg1', Nk, Nf, [])); % [Nk Nf Nt] 3D matrix, Aijk = Bij*Ck

Kernal = 2*pi*A_pump*(expterm .* repmat(kvect, [1 1 Nt])); % repmat projects kvect into t-space
% if memory is an issue, user can sacrifice readability to sandwich
% Kernal and expterm into the Integrand assignment. Shouldn't be necessary,
% since Kernal and expterm vanish after this function ends.

%Integrand=repmat(G, [1 1 Nt]) .* Kernal; % this reduces to G.*Kernal for Nt = 1.
%% GTH edit: September 4th, 2015: SushiTrain!
S = GetHankel_OffsetBeam_SushiTrain(kvect,xoffset,r_probe,1); % size(S) = [Nk, Nw, NFreq], where Nw = Nt.
S = permute(S,[1 3 2]); % Now size(S) = [Nk, NFreq, Nt], like the Kernal.

G = repmat(G, [1 1 Nt]); % Now size(G) = size(Kernal) = size(S).

% Integrand reduces to G.*Kernal.*S/(2*pi) for Nt = 1.
Integrand = G .* Kernal .* S; % divide by 2 pi, or not 2 pi?

% In DOUGHNUT, Integrand was G.*P.*S.*kvect, 
% where P.*kvect/(2*pi) = Kernal. I don't know where the (2*pi) was dropped
% when Joe went from TDTR_TEMP to TDTR_DOUGHNUT. Maybe S replaces it; that
% seems geometrically sensible.
end




        
