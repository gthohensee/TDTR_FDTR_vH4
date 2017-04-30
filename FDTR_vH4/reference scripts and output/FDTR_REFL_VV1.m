%Calculates the FDTR Reflectivity response
%In order to speed up the code, it is parallelized...the convention is...
%f is a vector of frequencies (Hz) to evaluate FDTR temperature response.
%kvectin (the fourier components of the in-plane temperature) are ROW vectors
function [deltaR,ratio,phi]=FDTR_REFL_VV1(f,lambda,C,t,eta,r_pump,r_probe,TCR,A_pump)

f = f(:)'; % frequencies Hz (convert to column vector if they aren't already)
kmax=1/sqrt(r_pump^2+r_probe^2)*10; %defines effective integration region

%Complex Temperature response as a function of frequency
[deltaT_model]=rombint_VV1(@(kvectin) FDTR_TEMP_VV1(kvectin,f,lambda,C,t,eta,r_pump,r_probe,A_pump),0,kmax,length(f)); 

%Reflectivity response (complex)
deltaR=TCR*deltaT_model;

%Ratio (X/Y)
ratio=real(deltaR)./imag(deltaR); 

%phase (in degrees)
phi=-angle(deltaR)*(180/pi); %phase delay, degrees

