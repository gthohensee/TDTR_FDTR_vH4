%% Define the directory and filenames
directory = '\\smb.oberon.ncf.wdc.com\EAMR\Gregory Hohensee\Data\TDTR'; % pwd is the current directory
yymmdd = '160629_Mindy_Rh_Y'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('\',yymmdd,'\'); % folder in directory where processed data is kept
savefolder = strcat('\',yymmdd,'_edit\'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

% get filenames in alphanumerical order
fnames = dir(strcat(datadir,'negboff*'));
nfiles = length(fnames); % number of data files
fstr = '';
for qq = 1:nfiles
    fstr = strcat(fstr, fnames(qq).name,'\n');
end
sprintf(fstr)

%% Process all FDTR files by subtracting the reference phase (leaked pump signal)
for k= 1:nfiles
  fname = fnames(k).name;
  
  data = dlmread(strcat(directory,datafolder,fname));
  
  xpos = data(:,1);
  Vin = data(:,3); Vout = data(:,4);
  Vin2 = Vin.^2;
  Vout2 = Vout.^2;
  AmpA = sqrt(Vin2+Vout2);
  PhiA = data(:,5);
  
  gfit(k,1:3) = fminsearch(@(X) SpotSize_FWHM_vH3(X,AmpA,xpos),[5 max(AmpA) 0]);
%     gfit =
%           2.81889709269829          116.421782255353        -0.484349515391256
%           2.70775499370759          108.950454810275        -0.363042519554749
%           2.70994626633191          114.551917301148        0.0880752211367002
    
    % From geometric mean: 2.56 um correlated spot.
end
gfit
%%
dlmwrite(strcat(savedir,'gaussian_fit_results.txt'), gfit, 'delimiter', '\t', 'precision', 4, 'newline', 'pc');
%% Any additional notes

% -------- END OF CODE --------
