%% Subfunction of TDTR_MAIN: Save and print figure from MAIN results
figsens = 202;
figfit = 203;

% Save the workspace
save(strcat(savedir,'Results_', fname(1:end-4),'.mat'))

% define the data and model arrays for the results figure
if ~doughnut
    switch sigfit
        case 1
            plot_data = Vin_data / Vin_data(Zind);
            plot_model = real(deltaR_model);
            plot_model = plot_model / (plot_model(Zind)/N); % see MANUALFIT for N
            ytext = 'normalized V(in)';
            fittext = 'FITVin_';
        case 2
            plot_data = Vout_data / Vout_data(Zind);
            plot_model = imag(deltaR_model);
            plot_model = plot_model / (plot_model(Zind)/N);
            ytext = 'normalized V(out)';
            fittext = 'FITVout_';
        case 3
            plot_data = Amp_data / Amp_data(Zind);
            plot_model = Amp_model;
            plot_model = plot_model / (plot_model(Zind)/N); % see MANUALFIT for N
            ytext = 'normalized amplitude';
            fittext = 'FITAmp_';
        case 4
            plot_data = phi_data;
            plot_model = phi_model;
            ytext = 'Phase (degrees)';
            fittext = 'FITPhase_';
        otherwise
            plot_data = ratio_data;
            plot_model = ratio_model;
            ytext = '-V(in)/V(out)';
            fittext = 'FIT_';
    end
else % beam offset
    switch sigfit
        case 1
            [~,I] = min(abs(offset));              % get index of r = 0 data point
            Norm = mean(abs(Vin_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
            plot_data = abs(Vin_data/Norm);      % remove negative sign, rescale to match DOUGHNUT.

            Vin_model = real(deltaR_model);        
            plot_model = abs(Vin_model) / max(abs(Vin_model));
            ytext = 'normalized V(in)';
            fittext = 'FITVin_boffset_'; 
        case 2
            [~,I] = min(abs(offset));              % get index of r = 0 data point
            Norm = mean(abs(Vout_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
            plot_data = abs(Vout_data/Norm);      % remove negative sign, rescale to match DOUGHNUT.

            Vout_model = imag(deltaR_model);        
            plot_model = abs(Vout_model) / max(abs(Vout_model));
            ytext = 'normalized V(out)';
            fittext = 'FITVout_boffset_'; 
        case 3
            plot_data = Amp_data ./ max(Amp_data);
            plot_model = Amp_model ./ max(Amp_model);
            ytext = 'Amplitude';
            fittext = 'FITAmp_boffset_';
        case 4
            plot_data = phi_data;
            plot_model = phi_model;
            ytext = 'Phase (degrees)';
            fittext = 'FITPhase_boffset_';
        otherwise
    end
    
end

% Create a figure for the last fit
figure(figfit)
clf;

if importdata
    w0 = sqrt(r_pump.^2 + r_probe.^2)/sqrt(2); % rescales x-axis on offset plot.
    if ~doughnut
        hd = semilogx(f_data/1e6,plot_data,'ko',...
                   'MarkerSize',8.0,'LineWidth',1);
        hold on;
        hm = semilogx(f_data/1e6,plot_model,'r-',...
                   'LineWidth',2);
    else
        hd = plot(offset*1e-6/w0 - V,plot_data / N,'ko',...
                   'MarkerSize',8.0,'LineWidth',1);
        hold on;
        hm = plot(offset*1e-6/w0,plot_model,'r-',...
                   'LineWidth',2);
    end

else % if importdata is false, then we didn't fit anything,
     % so there's no data points to plot -- just the model.
    if ~doughnut
        hm = loglog(f_data/1e6,plot_model,'r-',...
                   'LineWidth',2);
    else
        hm = plot(offset*1e-6/w0,plot_model,'r-',...
                   'LineWidth',2);
    end

end
fontsize = 16;

% define condtext
if P0 ~= -1
    condtext = sprintf('P = %0.1f GPa, Z = %0.2d',P0,Z);
elseif T0 ~= -1
    condtext = sprintf('T = %0.1f K, Z = %0.2d',T0,Z);
else
    condtext = sprintf('Z = %0.2d',Z);
end

% define axes
if doughnut
    switch sigfit
        case 1
            axis([min(offset/w0/1e6),max(offset/w0/1e6),-1,1.2]);
        case 2
            axis([min(offset/w0/1e6),max(offset/w0/1e6),-1,1.2]);
        case 3
            axis([min(offset/w0/1e6),max(offset/w0/1e6),-1,1.2]);
        case 4
            axis([min(offset/w0/1e6),max(offset/w0/1e6),-200,0]);
        otherwise
    end
else
    switch sigfit
        case 1
            set(gca, 'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 30 50]);
            axis([f_min/1e6 f_max/1e6 min(0.1,min(plot_data)) 2*max(plot_data)])
        case 2
            set(gca,'YScale','linear');
            set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
            axis([f_min/1e6 f_max/1e6 0 1.2]);
        case 3
            set(gca, 'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 30 50]);
            axis([f_min/1e6 f_max/1e6 min(0.1,min(plot_data)) 2*max(plot_data)])
        case 4
            axis([f_min/1e6 f_max/1e6 -90 0]);
        otherwise
    end
end


% define labels
if doughnut, xlabel('Offset x/w0', 'FontSize',fontsize); 
else xlabel('Frequency (MHz)','FontSize',fontsize); end

ylabel(ytext, 'FontSize',fontsize)
title(condtext,'FontSize',fontsize)

if ~doughnut
    set(gca, 'XTick', [0.1 0.2 0.5 1,2,5,10,20,50]);
    set(gca, 'XMinorTick', 'off');
    set(gca, 'YMinorTick', 'off');
end

set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',fontsize);
%% Provide a summary of final fit parameters in the figure.

% Warning: uses "ii" index from the analyze script's for loop!
% Compose fitstr and update LCTsol to match Xsol
XsolIJ = XguessIJ;
XsolIJ(:,1) = Xsol';
LCTEtext = genLCTEtext(LCTE, XsolIJ, stack(ii,:));

% write LCTEstr contents to a text box in the figure
pBox = annotation('textbox',[0.15,0.15,0.8,0.3]);
set(pBox,'Units','characters')
set(pBox,'HorizontalAlignment','left')
set(pBox,'FontSize',10)
set(pBox,'String',LCTEtext);
set(pBox,'LineStyle','none');

% also record the data file used
dBox = annotation('textbox',[0.15,0.85,0.8,0.07]);
set(dBox,'Units','characters')
set(dBox,'Interpreter','none')
set(dBox,'HorizontalAlignment','left')
set(dBox,'FontSize',10)
set(dBox,'String',{'Data file:'; fname(1:end)});
set(dBox,'LineStyle','none');

%% Save the figure to .fig and .eps files
saveas(figfit, strcat(savedir,fittext,fname(1:end-4),'.fig'))
print(figfit,'-depsc',strcat(savedir,fittext,fname(1:end-4),'.eps'))

if senseplot
    tag = input('Name the sensitivity plot: ','s');
    save(strcat(savedir,'/SENS_', tag, '.mat'))
    figure(figsens)
    saveas(figsens, strcat(savedir,'/SENS_',tag,'.fig'))
    print(figsens,'-depsc',strcat(savedir,'/SENS_',tag,'.eps'))
end