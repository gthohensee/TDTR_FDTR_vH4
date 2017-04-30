function [deltaR,ratio, phi]=FDTR_REFL_vH4(f, matparams,sysparams,A_pump,...
                                     intscheme,nnodes,offset)
%FDTR_REFL_vH - Calculates the Reflectivity Signal and Ratio at an array
%of discrete frequencies (for FDTR); no pulse train, no time delay.
%In order to speed up the code, it is parallelized...the convention is...
%frequencies "f" is COLUMN vector of desired frequencies for FDTR
%
% Syntax:  [deltaR,ratio]=FDTR_REFL_vH4(matparams,sysparams,A_pump,...
%                                     intscheme,nnodes,Fortran)
%
% Inputs (scalars unless told otherwise):  
%    f       - 1D VECTOR of modulation frequencies
%    A_pump  - pump intensity (W)...doesn't effect ratio
% SYSPARAMS %
%    r_pump  - pump 1/e2 radius (m)
%    r_probe - probe 1/e2 radius (m)
% MATPARAMS %
%    LCTE    - vertcat(lambda,C,t,eta);
%     lambda  - VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%     C       - VECTOR of specific heats (J/m3-K)
%     t       - VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%               you still have to enter something).  (m)
%     eta     - VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%               isotropic
%    BI      - TRUE if bidirectional heat flow is present.
%    n_toplayer - Specifies number of layers above the plane where
%                 heat is deposited in the bidirectional heat flow model.
%    TCR     - temperature coefficient of reflectivity  (1/K)...doesn't affect ratio
%    doughnut - TRUE if this is a beam offset calculation. Can combine with
%               BI for bidirectional beam offset.
% OTHER %
%    intscheme - Numerical integration algorithm.
%                0 = Legendre-Gauss, 1 = rombint, 2 = Simpson.
%    nnodes   - For Legendre-Gauss integration, number of nodes. Default
%               should be at least 35.
%    offset   - 1D vector of pump-probe relative offset distances (MICRONS)
%
% Outputs:
%    deltaR - Complex number VECTOR. Real part is model V(in), imaginary
%             part is model V(out).
%    ratio  - model -V(in)/V(out).
%    phi    - model phase angle between V(in) and V(out).
%
% Other m-files required: FDTR_TEMP_vH4.m
% Subfunctions: none
% MAT-files required: none
%
% See also: The TDTR_V4 package.

% Author: Gregory Hohensee / built from Joseph Feser's TDTR_REFL_V4.
% U. of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history: 12-Sep-2012 - J. Feser's REFL_V4 published
%                   2013 - added bidirectional heat flow, REFL_V4B.
%                   25-Mar-2014 - standardized comments into TDTR_vH package.
%                   26-Mar-2014 - use LCTE instead of lambda,C,t,eta. 
%                                 Allow various integration schemes.
%                   8-Apr-2014 - harmonized with TTM version
%                   14-July-2014 - vH2. No change.
%                   13-Nov-2014 - Merged with beam offset code.
%                   21-Jan-2015 - accepts offset as vector.
%                   28-Jan-2015 - error messages now block rombint and
%                                 Simpson int from choking on a variable
%                                 pump spot size. Only works with L-G int.
%                   17-Feb-2015 - vH3.
%                   1-Sep-2015  - FDTR, no spot size variation.
%                   11-Jul-2016 - vH4.
%------------- BEGIN CODE --------------
%% check input
if nargin < 6
    offset = 0; 
    if doughnut, warning('Assuming zero offset in REFL_vH4.'); end 
end
offset = offset*1e-6; % correct units for TEMP_DOUGHNUT

if nargin < 5, nnodes = 35; end
if nargin < 4, intscheme = 0; end
if nargin < 3, A_pump = 10e-3; end % arbitrary positive value
if nargin < 2, error('Insufficient parameters for FDTR_REFL_vH4.'); end

%% unpack input
PL = [0 length(sysparams) 0 length(matparams) 0];

% pump-probe system parameters
%sysparams = {r_pump r_probe};
if PL(2) < 1, error('pump spot size not specified.'); else r_pump = sysparams{1}; end
if PL(2) < 2, r_probe = r_pump;
    warning('Defaulting to r_probe = r_pump.'); else r_probe = sysparams{2}; end

% material parameters
%matparams = {LCTE aniso BI n_toplayer TCR};
if PL(4) < 1, error('LCTE not specified'); else LCTE = matparams{1}; end
%if pl(4) < 2, aniso = zeros(size(LCTE(4,:))); else aniso = matparams{2}; end
if PL(4) < 3, BI = 0; else BI = matparams{3}; end
if PL(4) < 4, if BI, error('n_toplayer not specified'); else n_toplayer = 0; end 
   else n_toplayer = matparams{4}; end
if PL(4) < 5, TCR = 1e-4; else TCR = matparams{5}; end
if PL(4) < 6, doughnut = 0; else doughnut = matparams{6}; end

%% MATLAB port, bi-directional anisotropic (!) thermal model.
% That is, I see no reason why this can't also handle anisotropic cases.

    if intscheme == 2 % Simpson integration
        % implement identically to RW scripts, which use Simpson with
        % nnodes = 35.
        cap2 = 2;
    else % go with JPF / prior factors.
        cap2 = 1.5;
    end
    
    kmin = 0;
    kmax=cap2/sqrt(min(r_pump)^2+r_probe^2);
    % RW uses cap2 = 2 instead of 1.5 for his REFL/TEMP scripts,
    % which use Simpson integration exclusively. 2x catches higher k-modes,
    % is more accurate assuming comparable node density over k-space.
    % For TTM, k(RW) = 4*pi^2*k(JPF), where k(JPF) is in use in TEMP_V4.
    
    if BI
        dT1u = zeros(1,length(f))';
        dT1d = zeros(1,length(f))';
    end
    dT1=zeros(1,length(f))';
    
% loop through pump-probe offset positions, because FDTR_TEMP_DOUGHNUT and
% associated Hankel calculations aren't parallelized in offset. Not sure
% how to do that, and we're only looping ~10 elements in offset anyway.
deltaR = zeros(length(f),length(offset));
ratio = deltaR;
phi = deltaR;

for jj = 1:length(offset)
    xo = offset(jj); % still fine if offset = 0 (time-delay TDTR)
    
    if BI
        totlayers = length(LCTE(1,:));

        % Split the material parameters into "up" and "down" sections,
        % relative to the heater and thermometer at interface of n_toplayer
        % and n_toplayer+1.
        for j = n_toplayer:-1:1
            LCTEu(:,n_toplayer+1-j) = LCTE(:,j);
        end
        LCTEd = LCTE(:,n_toplayer+1:totlayers);
        
        switch intscheme
            case 0 % Legendre-Gauss, allows variable pump spot size.
                [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
                % calculate the heat flow up, layer n_toplayer up to layer 1
                if doughnut
                    I1u = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump,xo);
                else
                    I1u = FDTR_TEMP_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump);
                end
                dT1u = weights'*I1u; % the k-space integral
                
                % calculate the heat flow down, layer n_toplayer+1 down to bottom
                if doughnut
                    I1d = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump,xo);
                else
                    I1d = FDTR_TEMP_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump);
                end
                dT1d = weights'*I1d; % the k-space integral
                
            case 1 % rombint [had a +/-f typo here before Nov. 13, 2014.]
                if doughnut
                    dT1u=rombint_VV3(@(kvect) FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump,xo),0,kmax,1);
                    dT1d=rombint_VV3(@(kvect) FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump,xo),0,kmax,1);
                else
                    dT1u=rombint_VV3(@(kvect) FDTR_TEMP_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump),0,kmax,1);
                    dT1d=rombint_VV3(@(kvect) FDTR_TEMP_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump),0,kmax,1);
                end
            case 2 % Simpson integration
                kvect = linspace(0,kmax,nnodes)';
                if doughnut
                    I1u = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump,xo);
                    I1d = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump,xo);
                else
                    I1u = FDTR_TEMP_vH4(kvect,f,LCTEu,r_pump,r_probe,A_pump);
                    I1d = FDTR_TEMP_vH4(kvect,f,LCTEd,r_pump,r_probe,A_pump);
                end
                dT1u=SimpsonInt(kvect,I1u);
                dT1d=SimpsonInt(kvect,I1d);
            otherwise
                error('integration scheme not properly specified. See analyze template.');
        end
        % make the parallel sum of the temperatures for up and down
        dT1 = dT1u .* dT1d ./ (dT1u + dT1d);
        
    else % BI = 0, not bidirectional.    
        switch intscheme
            case 0 % Legendre-Gauss [compatible with time delay dependent spot size]
                [kvect,weights]=lgwt_V4(nnodes,kmin,kmax); %computes weights and node locations...
                
                if doughnut % TRUE if this is a beam offset measurement
                    I1 = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump,xo);
                else
                    I1 = FDTR_TEMP_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump);
                end
                
                dT1 = weights'*I1;
            case 1 % rombint
                if doughnut
                    
                    dT1=rombint_VV3(@(kvect) FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump,xo),0,kmax,1);
                else
                    dT1=rombint_VV3(@(kvect) FDTR_TEMP_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump),0,kmax,1);
                end
                
            case 2 % Simpson
                kvect = linspace(0,kmax,nnodes)';
                if doughnut
                    I1 = FDTR_TEMP_DOUGHNUT_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump,xo);
                else
                    I1 = FDTR_TEMP_vH4(kvect,f,LCTE,r_pump,r_probe,A_pump);
                end
                dT1=SimpsonInt(kvect,I1);
                
            otherwise
                error('integration scheme not properly specified. See analyze template.');
        end
    end
    
    % dT1 should be a 1D column vector in frequency space up to this line.
    dT1 = dT1(:)'; % switch to row vector
    
    % columns are offset (the jj index), rows are frequency space.
    %Reflectivity response (complex)
    deltaR(:,jj)=TCR*dT1;
    
    %Ratio (X/Y)
    ratio(:,jj)= real(deltaR(:,jj))./imag(deltaR(:,jj)); 

    %phase (in degrees)
    phi(:,jj) = -angle(deltaR(:,jj))*(180/pi); %phase delay, degrees

end % end for loop for offset positions
end