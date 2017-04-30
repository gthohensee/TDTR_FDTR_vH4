%% Define the directory and filenames
directory = 'pwd'; % pwd is the current directory
yymmdd = 'yymmdd'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('\',yymmdd,'\'); % folder in directory where processed data is kept
savefolder = strcat('\',yymmdd,'_edit\'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

% get filenames in alphanumerical order
fnames = dir(strcat(datadir,'boffset*'));
nfiles = length(fnames); % number of data files
fstr = '';
for qq = 1:nfiles
    fstr = strcat(fstr, fnames(qq).name,'\n');
end
sprintf(fstr)
clear gfit

%% Process all FDTR files
for k= 1:nfiles
  fname = fnames(k).name;
  
  data = dlmread(strcat(directory,datafolder,fname));
  
  xpos = data(:,1);
  Vin = data(:,3); Vout = data(:,4);
  Vin2 = Vin.^2;
  Vout2 = Vout.^2;
  AmpA = sqrt(Vin2+Vout2) - .2;
  PhiA = data(:,5);
  
  gfit(k,1:3) = fminsearch(@(X) SpotSize_FWHM_vH3(X,AmpA,xpos),[2 max(AmpA) 0]);
  
    % 4/27/2016. Assume 3 microns isotropic 532 nm spot size (10 mW laser),
    % from knife edge measurements.
    % 4/29/2016: 20x is basically 1.98, aka 2.0 micron correlated spot size.
    % 
end

% pull out pump spot size from correlated spot size
gfit(:,5) = abs(sqrt(2*gfit(:,1).^2 - 1.5^2)); % assuming probe spot size of 3/2 microns (FDTR 532 nm) at 20x

% print gfit
gfit

% write gfit to file
dlmwrite(strcat(datadir,'gaussian_fit_results.txt'), gfit, 'delimiter', '\t', 'precision', 4, 'newline', 'pc');

%% Any additional notes

% -------- END OF CODE --------
