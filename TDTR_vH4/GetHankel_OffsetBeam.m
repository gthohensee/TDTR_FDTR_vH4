function S = GetHankel_OffsetBeam(kvect,xo,ws,A)

%kvect is Nk x NFreq matrix
k=kvect(:,1); % k is Nk column vector

% Add this for wpump varying across tdelay.
% Not guaranteed to work if wpump != wprobe?
ws = ws(:)'; % forces ws to be a row vector.

kbar = k*ws/sqrt(2); % now kbar is Nk x N(ws) matrix.
xbar = sqrt(2)*xo/ws; % xo is scalar, xbar now has dimension wpump.

if nargin ==3
    A =1;
end

prefactor = A/pi*exp(-(xbar^2+pi^2*kbar.^2));
x2 = xbar^2;
p = OffcenterBeam;
Nmax = length(p(:,1))-1;
Nk = length(kbar);

summand = zeros(size(k));
sigma = zeros(size(k));
for n = 0:Nmax
    PP=polyval(p(n+1,:),kbar);
    %PP=PP.*(~isinf(PP)); %protects against numerical instability:  DANGER
    summand = x2^n/factorial(n)^2*PP;
    sigma = sigma + summand;
    
%     %for debugging only
%     S = prefactor.*sigma;
%     plot(k,S)
%     hold on
%     pause(1)


end

S = (prefactor.*sigma)*ones(1,length(kvect(1,:)));

    