%% Define the directory and filenames
directory = 'pwd'; % pwd is the current directory
yymmdd = 'yymmdd'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('\',yymmdd,'\'); % folder in directory where processed data is kept
savefolder = strcat('\',yymmdd,'_edit\'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

% get filenames in alphanumerical order
fnames = dir(strcat(datadir,'*_1.txt'));    % pump-probe FDTR data
refnames = dir(strcat(datadir,'*ref.txt')); % reference pump phase data (block probe, leak pump)

nfiles = length(fnames); % number of data files

fstr = ''; refstr = '';
for qq = 1:nfiles
    fstr = strcat(fstr, fnames(qq).name,'\n');
    refstr = strcat(refstr,refnames(qq).name,'\n');
end
sprintf(fstr)
sprintf(refstr)


%% Process all FDTR files by subtracting the reference phase (leaked pump signal)
for k=1:nfiles
  fname = fnames(k).name;
  refname = refnames(k).name;
  
  data = dlmread(strcat(directory,datafolder,fname));
  refdata = dlmread(strcat(directory,datafolder,refname));
  
  AmpA = data(:,4);
  PhiA = data(:,5);
  refPhi = refdata(:,5);

  deltaR_corrected = AmpA .* exp(sqrt(-1) * (pi/180) * (PhiA - refPhi - 180)); 
  % why is it not just PhiA - refPhi?
  % Because I'm an idiot: the TCR thermoreflectance coefficient is NEGATIVE
  % for Au, but the default value I've been using is POSITIVE (for Al). So
  % that flips the sign of deltaR, hence also adds 180 degree phase shift.
  
  cVin = real(deltaR_corrected);
  cVout = imag(deltaR_corrected);
  cPhi = (180/pi) * angle(deltaR_corrected); 

  % cdata is corrected data
  clear cdata
  cdata(:,1) = data(:,1); 
  cdata(:,2) = cVin;
  cdata(:,3) = cVout;
  cdata(:,4) = data(:,4);
  cdata(:,5) = cPhi;
  
  % save to file
  fileprefix = fname(1:end-4);
  file_shifted=strcat(savedir,fileprefix,'_shifted.txt');
  
  dlmwrite(file_shifted, cdata, 'delimiter', '\t', 'precision', 4, 'newline', 'pc')
  fprintf('exporting file to:\n');
  fprintf('%s\n',file_shifted);
end
%%

%% Any additional notes

% -------- END OF CODE --------
