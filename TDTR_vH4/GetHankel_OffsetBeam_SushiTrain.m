function S = GetHankel_OffsetBeam_SushiTrain(kvect,xo,ws,A)
% SushiTrain is for use in TDTR_TEMP_SUSHITRAIN, spinning thermal model.
% This was a test which was never validated; adapt with caution.

%kvect is Nk x NFreq matrix
Nk = length(kvect(:,1));
NFreq = length(kvect(1,:));
Nx = length(xo);

k=kvect(:,1); % k is Nk column vector

xbar = sqrt(2)*xo(:)'/ws; % xo is array, xbar is now a row vector
xbar = repmat(xbar,Nk,1); % now xbar is Nk x N(ws), with identical entries along Nk.

kbar = k*ws/sqrt(2);
kbar = repmat(kbar,1,Nx);

if nargin == 4
    A =1;
end

x2 = xbar.^2;
prefactor = A/pi*exp(-(x2+pi^2*kbar.^2));

p = OffcenterBeam;
Nmax = length(p(:,1))-1;

sigma = zeros(size(kbar));
for n = 0:Nmax
    PP=polyval(p(n+1,:),kbar); % PP is a matrix of size kbar
    %PP=PP.*(~isinf(PP)); %protects against numerical instability:  DANGER
    summand = x2.^n/factorial(n)^2 .* PP; % x2.^n and PP have same size.
    sigma = sigma + summand;
    
%     %for debugging only
%     S = prefactor.*sigma;
%     plot(k,S)
%     hold on
%     pause(1)


end

S1 = prefactor.*sigma; % has dimensions (Nk, Nw), where Nw = Nt is the number of wpump (tdelay) elements.
S = repmat(S1, [1 1 NFreq]); % now S has dimension (Nk, Nx, NFreq)

    