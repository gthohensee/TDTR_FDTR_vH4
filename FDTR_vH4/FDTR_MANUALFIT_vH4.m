function [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj,N,V]=...
    FDTR_MANUALFIT_vH4(XguessIJ,datparams,...
                       sysparams, calparams, matparams, Tparams)
%FDTR_MANUALFIT_vH4 - Manually fit thermal model to FDTR ratio data.
%This program lets you iteratively fit the thermal model to the ratio data
%by varying the LCT thermal parameters specified in XguessIJ. It will tell
%you the residual deviation (crude goodness-of-fit) Z, allow you to
%generate a sensitivity plot at any time, and do thermal modeling with
%anisotropic unidirectional or bidirectional heat flow. It also handles
%temperature dependence in a self-consistent manner with steady-state and
%per-pulse heating through calls to TDTR_TDEP_vH4.m.
%
%It will NOT allow you to:
% ** reduce the transducer thickness below the absorption length if you
%    model an absorption layer, 
% ** compute the goodness-of-fit weighted by the sensitivities. 
%    It'd take too long to generate sensitivities every time you change 
%    the fit.
%
% Syntax:  [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj]=...
%    TDTR_MANUALFIT_vH4(XguessIJ,datparams,...
%                       sysparams, calparams, matparams, Tparams)
%
% Inputs:
%    XguessIJ  - Mx3 matrix: each row represents a fit parameter Xguess,
%                so [Xguess i j]. (i,j) designates Xguess => LCTE(i,j).
%    datparams - {tdelay ratio_data datadir offset}
%    sysparams - {tau_rep f r_pump r_probe}
%    calparams - {Zind sigfit intscheme nnodes consider_error LCTE_err T0_err}
%    matparams - {LCTE aniso BI n_toplayer TCR doughnut}
%    Tparams   - {T0, T_LCTE, A_pump, A_probe, absC, perpulse, jabs, jtrans}
%[See INITIALIZE_CELLPARAMS_vH4.m for details on the params inputs.]
%
% Outputs:
%    Xsol        - Final values for the fitted parameters.
%    Z           - Goodness-of-fit. In Joe's words: "Typically, this is 
%                  the sum of the squares of the residuals, but you might 
%                  want to weight these by the sensitivity, particularly 
%                  if you don't intend to calculate the errorbars!"
%    deltaR      - Complex number array. Real part is the model V(in),
%                  imaginary part is the model V(out). Represents change
%                  in reflectance from pump heating.
%    ratio_model - Ratio signal -V(in)/V(out) from the thermal model.
%    LCTE        - The values in Xsol may be tied to the anisotropy eta,
%                  or they could affect T_adj, which affects values in
%                  LCTE generally. So, output a new LCTE.
%    T_adj       - T_adj = T0 + dTss + dTpp, the "actual" temperature
%                  in Kelvin adjusted for steady-state (and per-pulse
%                  if perpulse is TRUE) heating.
%
% Example:
%    --
%
% Other m-files required: FDTR_REFL_vH4.m, FDTR_TEMP_vH4.m, 
%                         TDTR_TDEP_vH4.m, SS_Heating_vH4.m
% Subfunctions: none
% MAT-files required: none
%
% See also: TDTR_FIT_TTM_vH4.m, TDTR_MANUALFIT_TTM_vH4.m, TDTR_FIT_vH4.m

% Author: Gregory Hohensee
% Acknowledgement: built from TDTR_FIT_V4B, my bi-directional tweak to the
% original TDTR_FIT_V4.m, written by the great and powerful Joseph P.
% Feser.
% University of Illinois at Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% Revision history:  25-Mar-2014 - written
%                    8-Apr-2014 - more comments, harmonized with TTM
%                    14-July-2014 - vH2. No changes.
%                    13-Nov-2014 - now allows beam offset (doughnut),
%                                  only for V(out) and no sensitivity
%                                  calculations.
%                    27-Jan-2014 - examined, no code changes.
%                    17-Feb-2015 - vH3, now calls "senseplot_offset_vH4"
%                                  instead of "senseplot_vH2" for
%                                  sensitivities.
%                    1-Sep-2015  - FDTR.
%                    11-Jul-2016 - vH4.
%------------- BEGIN CODE --------------
%% Check input parameters, assign defaults, errors, warnings as necessary
FDTR_INITIALIZE_CELLPARAMS_vH4; % unpacks and re-packs cellparams

% TEMPORARY, FOR 151216/7 ONLY: subtract 0.05-0.1 uV from Amplitude
% model to account for background noise in boffset (no chopper).
% Amp_data = Amp_data - 0.1;

%% Assign fit variables X according to their XguessIJ index.
% XguessIJ = [X,i,j] indexes X across LCTE(i,j)
Xij = XguessIJ(:,2:3); % for easier comparison to TDTR_FIT_vH.m
nrows = length(Xij(:,1));

% iterate over fit parameters
for x = 1:nrows, LCTE(Xij(x,1),Xij(x,2)) = XguessIJ(x,1); end 
%% User input loop
done = 0; 
N = 1; % N is for adjusting the normalization factor when fitting by V(in) or V(out)
V = 0; % V is for adjusting the lateral position of the offset data.
w0 = sqrt(r_pump.^2 + r_probe.^2)/sqrt(2); % rescales x-axis on offset plot.

while done ~= 1
    %% Self-consistent steady-state (and optionally per-pulse) heating
    % NOT UPDATED FOR FDTR; this is just leftover from TDTR_vH3.
    if T0 ~= -1
        error('I did not update self-consistent T adjustment to the FDTR setup. Go back to T0 = -1.')
        [dTss, dTpp, LCTE] = TDTR_TDEP_vH4(matparams,sysparams,...
                                               Tparams,intscheme,nnodes);
        matparams{1} = LCTE; % update                                   
        fprintf('T0 = %0.2f K, dTss = %0.2f, dTpp = %0.2f\n',T0,dTss,dTpp)
    end
    
    %% Compute model at zero beam offset (conventional TDTR)
    [deltaR_model,ratio_model,phi_model]=FDTR_REFL_vH4(f_data,matparams,sysparams,...
                                                A_pump,intscheme,nnodes,0);
    Vin_model = real(deltaR_model);
    Vout_model = imag(deltaR_model);
    Amp_model = sqrt(Vout_model.^2 + Vin_model.^2);
    max(Amp_model)
    if ~doughnut
        switch sigfit
            case 1 % V(in) fit
                % Construct normalized V(in) model and data,
                % relative to its value at ZF frequency.

                Vin_model_ZF = Vin_model(Zind) / N;
                nVin_model = Vin_model / Vin_model_ZF;

                Vin_data_ZF = Vin_data(Zind);
                nVin_data = Vin_data / Vin_data_ZF;
            case 2 % V(out) fit
                Vout_model_ZF = Vout_model(Zind) / N;
                nVout_model = Vout_model / Vout_model_ZF;

                Vout_data_ZF = Vout_data(Zind);
                nVout_data = Vout_data / Vout_data_ZF;
            case 3
                Amp_model_ZF = Amp_model(Zind) / N;
                nAmp_model = Amp_model / Amp_model_ZF;

                Amp_data_ZF = Amp_data(Zind);
                nAmp_data = Amp_data / Amp_data_ZF;
            case 4
                %phi_model is phase delay, degrees, computed by REFL
            otherwise % ratio fit
                % ratio_model computed by REFL
        end
    else
    % Do beam offset thermal model.
    % input: offset, Nx, REFL params. output: deltaR and ratio vs xvect
        Y0 = imag(deltaR_model);  % this from a REFL call with xoffset = 0.
        X0 = real(deltaR_model);
        
        
        % fyi -- normally no need to model past 3.9x spot size.
        %Nx = 10;
        %xvect = linspace(0,max(offset),Nx)'; % need to model just one side of the data, by symmetry
        xvect = offset; % brute solution to calculating residual
        
        % calculate model from each offset point in xvect.
        [deltaR_modelx, ratio_modelx, phi_modelx]=FDTR_REFL_vH4(f_data,matparams,sysparams,...
                                                   A_pump,intscheme,nnodes,xvect);
	    % deltaR and ratio may be a matrix (frequencies, offset). Pull out
        % first (only) frequency element, transpose into column. Beam
        % offset is taken at a single frequency.
        Sig_Amp = deltaR_modelx(1,:)'; 
        Sig_ratio = ratio_modelx(1,:)';
        Sig_phi = phi_modelx(1,:)';
        
        %reconstruct the rest by symmetry
        %xvect = [-flipud(xvect(2:end));xvect];
        %Sig_Amp = [flipud(Sig_Amp(2:end));Sig_Amp];
        %Sig_ratio = [flipud(Sig_ratio(2:end));Sig_ratio];
        %Sig_phi = [flipud(Sig_phi(2:end));Sig_phi];
        
        Vout_model = interp1(xvect,imag(Sig_Amp),offset,'spline','extrap');
        Vin_model = interp1(xvect,real(Sig_Amp),offset,'spline','extrap');
        deltaR_model = Vin_model + 1i*Vout_model;
        ratio_model = -Vin_model ./ Vout_model;
        Amp_model = Sig_Amp;
        
        [~,I] = min(abs(offset));
        
        
        
        switch sigfit
            case 1 % Vin
                NormX = mean(abs(Vin_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
                nVin_data = abs(Vin_data/NormX);      % remove negative sign, rescale to match DOUGHNUT.
            case 2 % Vout
                NormY = mean(abs(Vout_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
                nVout_data = abs(Vout_data/NormY);      % remove negative sign, rescale to match DOUGHNUT.
            case 3 % amplitude
                NormA = mean(abs(Amp_data(I-1:I+1)));
                nAmp_data = abs(Amp_data/NormA);
            case 4 % phase
                phase_offset = Sig_phi(I) - N*phi_data(I);
                phi_model = angle(exp(sqrt(-1)*(Sig_phi - phase_offset)*pi/180))*180/pi;
                fprintf('Note: automatically applied phase offset of %0.2f degrees to the model to fit x = 0',phase_offset);
            otherwise % ratio
        end
        
        nVout_model = abs(Vout_model/Y0);
        nVin_model = abs(Vin_model/X0);
        nAmp_model = Amp_model / max(Sig_Amp);
    end
    %% Update the data and fit comparison figure
    figure(10)
     
    switch sigfit
        case 1
            loglog(f_data/1e6,nVin_data,'ob',f_data/1e6,nVin_model,'r');
            
            ylabel('normalized V(in)','FontSize',20);
            axis([min(f_data/1e6) 2*max(f_data/1e6) min([nVin_data;nVin_model])/2 1.4])
        case 2
            if ~doughnut
                semilogx(f_data/1e6,nVout_data,'ob',f_data/1e6,nVout_model,'r');

                ylabel('normalized V(out)','FontSize',20);
                axis([1e-10 10e-9 min([0;nVout_data;nVout_model]) 1.4])
                set(gca,'YMinorTick','on')
                set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
            else % beam offset
                axis([-4 4 -0.2 1.4]);
                plot(offset*1e-6/w0 - V, nVout_data/N,'ok',...
                     xvect*1e-6/w0,      nVout_model,'r-','LineWidth',2);
                xlabel('x/w_0','FontSize',20)
                ylabel('Offset normalized V(out)','FontSize',20)
            end
        case 3 % amplitude
            if ~doughnut
                loglog(f_data/1e6,nAmp_data,'ob',f_data/1e6,nAmp_model,'r');
                ylabel('normalized amplitude','FontSize',20);
                axis([min(f_data/1e6) 2*max(f_data/1e6) min([nAmp_data;nAmp_model])/2 1.4])
            else
                axis([-4 4 -0.2 1.4]);
                plot(offset - V, nAmp_data/N,'ok',...
                     xvect,      nAmp_model, 'r-','LineWidth',2);
                xlabel('x/w_0','FontSize',20)
                ylabel('Offset normalized amplitude','FontSize',20)
            end
        case 4 % phase
            if ~doughnut
                semilogx(f_data/1e6,phi_data,'ob',f_data/1e6,phi_model,'r');
                ylabel('phase (degrees)','FontSize',20);
                axis([min(f_data/1e6) 2*max(f_data/1e6) -90 0])
            else
                axis([-4 4 -90 0]);
                plot(offset*1e-6/w0 - V, phi_data, 'ok',...
                     xvect*1e-6/w0,      phi_model,'r-','LineWidth',2);
                xlabel('x/w_0','FontSize',20)
                ylabel('Offset phase (degrees)','FontSize',20)
            end
        otherwise
            loglog(f_data/1e6,ratio_data,'ob',f_data/1e6,ratio_model,'r');
            ylabel('Ratio','FontSize',20); 
            axis([min(f_data/1e6) 2*max(f_data/1e6) min([ratio_data;ratio_model])/2 max([ratio_data;ratio_model])*2])
            
    end
    if sigfit ~= 4 && sigfit ~= 3
        set(gca, 'YTick',      [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca, 'YTickLabel', [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca,'YMinorTick','off')
    end
    
    
    set(gca,'FontSize',20);
    set(gca, 'TickLength' , [.02 .02]);
    if ~doughnut
        xlabel('Frequency (MHz)','FontSize',20);
        set(gca, 'XTick',      [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca, 'XTickLabel', [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]);
        set(gca, 'XMinorTick','off')
    end
    
    
    % goodness-of-fit Z from fractional residuals
    switch sigfit
        % For V(in) and V(out) fits, it may be wiser to weight the 
        % residuals by the amount of time delay away from ZF, 
        % where V(in) is pinned to unity.
        case 1
            res=(1-(nVin_model(Zind:length(nVin_model)) ...
                    ./nVin_data(Zind:length(nVin_data))) ).^2;
        case 2
            if doughnut
                res = (nVout_model - nVout_data).^2; % fractional error may overvalue the tails?
            else
                res=(1-(nVout_model(Zind:length(nVout_model)) ...
                        ./nVout_data(Zind:length(nVout_data))) ).^2;
            end
        case 3
            if doughnut
                res = (nAmp_model - nAmp_data).^2;
            else
                res=(1-(nAmp_model(Zind:length(nAmp_model)) ...
                        ./nAmp_data(Zind:length(nAmp_data))) ).^2;
            end
        case 4
            if doughnut
                res = ((phi_model - phi_data)./(phi_data)).^2;
            else
                res=(1-(phi_model(Zind:length(phi_model)) ...
                        ./phi_data(Zind:length(phi_data))) ).^2;
            end
        otherwise % ratiofit
        res=(1-(ratio_model(Zind:length(ratio_model)) ...
                ./ratio_data(Zind:length(ratio_data))) ).^2;
    end 
    Z = sum(res)
    
    %% Inform user of current fit parameters
    
    % Define symbols and units for all possible fit parameters
    tag = {'L','C','t','eta'};
    units = {'W/m-K', 'J/cm^3-K', 'nm', '(Lx/Lz)'};
    scale = [1 1e-6 1e9 1];
    
    % Report current fit parameters
    Mat = sprintf('Current material fit parameters LCTE(i,j):');
    for x = 1:nrows % iterate through fit parameters
        itag = Xij(x,1);
        Mat = char(Mat,sprintf('%s(%i) = %0.4f %s',...
                   tag{itag},Xij(x,2),scale(itag)*LCTE(Xij(x,1),Xij(x,2)),units{itag}));
    end
    Mat
    
    %% get and clean input
    done = input('Enter 1 if done, 2 for a sensitivity plot, 3 to rescale; else hit "Enter": ');
    if isempty(done) 
        done = 0;
    else
        if sum(done == [0 1 2 3]) == 0
            done = 0; fprintf('Hey! Invalid input. Go home, you are drunk.\n')
        end
    end
    %% execute user choice
    switch done
        case 3 % tweak the normalization factor
            switch sigfit
                case 1
                    fprintf('Pick the normalization point on the plot...\n')
                    [V,N] = ginput(1);
                    fprintf('You picked [%0.2f %0.2f]. Calculating new fit...\n',V,N);
                case 2
                    fprintf('Pick the normalization point on the plot...\n')
                    [V,N] = ginput(1);
                    fprintf('You picked [%0.2f %0.2f]. Calculating new fit...\n',V,N);
                case 3
                    fprintf('Pick the normalization point on the plot...\n')
                    [V,N] = ginput(1);
                    fprintf('You picked [%0.2f %0.2f]. Calculating new fit...\n',V,N);
                case 4 % phase
                    fprintf('Pick the normalization point on the plot...\n')
                    [V,N] = ginput(1);
                    N = N ./ phi_data(I);
                    fprintf('You picked [%0.2f %0.2f xPhi0]. Calculating new fit...\n',V,N);
                otherwise
                    fprintf('You cannot normalize the ratio!\n')
            end
        case 2 % sense plot
            %if doughnut 
            %    warning('Sorry, I still need to code the beam offset sensitivity calculation for non-Vout.');
            %    break;
            %end
            
            [S_LCTE, S_sys] = FDTR_senseplot_offset_vH4(datparams,sysparams, calparams, matparams, Tparams);
            savesens = input('Enter 1 to save sensitivity plot to datadir: ');
            if savesens == 1
                tag = input('Name the sensitivity plot: ','s');
                save(strcat(datadir,'/SENS_', tag, '.mat'))
                fignum=202;
                figure(fignum)
                saveas(fignum, strcat(datadir,'/SENS_',tag,'.fig'))
                print(fignum,'-depsc',strcat(datadir,'/SENS_',tag,'.eps'))
            end
        case 0 % manual fitting
            for x = 1:nrows % iterate through fit parameters
                itag = Xij(x,1);
                temp = input(sprintf('Adjust parameter %s(%i): ',...
                                        tag{itag},Xij(x,2)));
                if ~isempty(temp)
                    XguessIJ(x,1) = temp / scale(itag);
                    
                    % if adjusting conductivity, then eta may
                    % change, depending on aniso.
                    if Xij(x,1) == 1 && aniso(Xij(x,2))
                        LCTE(4,Xij(x,2)) = LCTE(4,Xij(x,2)) * LCTE(1,Xij(x,2)) / temp;
                    end

                    % having used the old LCTE(1,j) to update eta if necessary,
                    % can now update LCTE(i,j).
                    LCTE(Xij(x,1),Xij(x,2)) = temp / scale(itag);
                end
                
            end % end iteration through fit parameters
    end % end switch for execution of user's decision
    
    matparams{1} = LCTE; % update matparams.
end % end user input loop

% Assign final temperature adjusted for SS and PP heating.
if T0 ~= -1, T_adj = T0 + dTss + dTpp; else T_adj = T0; end

% Assign Xsol values
Xsol = XguessIJ(:,1);
end % end program
