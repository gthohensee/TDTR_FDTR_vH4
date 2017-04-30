function [Integrand,G]=FDTR_TEMP_vH4(kvectin,freq,LCTE,r_pump,r_probe,A_pump)
%FDTR_TEMP_vH3 - Computes frequency domain average temperature response to 
%periodic gaussian pump beam, probed by another gaussian beam.
%
% This program is vectorized/optimized to handle all frequencies 
% simultaneously (f is a ROW vector)
%
% Syntax:  [Integrand,G]=FDTR_TEMP_vH4(kvectin,freq,LCTE,r_pump,r_probe,A_pump)
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
% See also: FDTR_REFL_vH4.m

% Author: Joseph P. Feser
% --
% email: --
% Website: TDTR_V4 package: http://users.mrl.illinois.edu/cahill/tcdata/tcdata.html
% September 2012; 
% Revision history: 12-September-2012 - published as TDTR_TEMP_V4.m
%                   25-Mar-2014 - [header comments reformatted by Greg Hohensee]
%                   8-Apr-2014 - harmonized with TTM version
%                   14-July-2014 - vH2. No change.
%                   17-Feb-2015  - vH3. No change.
%                   1-Sep-2015   - converted from TDTR to FDTR. No variable
%                                  pump spot size.
%                   11-Jul-2016  - vH4. No change.
%------------- BEGIN CODE --------------
%% variables from previous versions of this code.
% Unpacked from my other data structures so I don't need to edit the code.
lambda = LCTE(1,:);
C = LCTE(2,:);
t = LCTE(3,:);
eta = LCTE(4,:);

Nfreq=length(freq);
kvect=kvectin(:)*ones(1,Nfreq); % Nfreq rows of (0,...,kmax) column vectors
Nlayers=length(lambda); %# of layers

Nint=length(kvectin); %# of different spatial frequencies to calculate for

%% Joe's code
%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

ii=sqrt(-1);
alpha=lambda./C;
omega=2*pi*freq(:)';
q2=ones(Nint,1)*(ii*omega./alpha(Nlayers));
kterm2=4*pi^2*kvect.^2; % GTH: is this term why kmax ~ 1.5 instead of 4*pi when defining kvectin?

% Looks like the n = Nlayers start point of the calculation.
un=sqrt(eta(Nlayers)*kterm2+q2);
gamman=lambda(Nlayers)*un;
Bplus=zeros(Nint,Nfreq);
Bminus=ones(Nint,Nfreq);

if Nlayers~=1
    for n=Nlayers:-1:2
        % n only appears in this loop as (n-1). Count is Nlayers to 1.
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

Kernal=2*pi*A_pump*exp(-pi^2*(r_pump^2+r_probe^2)/2*kvect.^2).*kvect; %The rest of the integrand
Integrand=G.*Kernal;

end