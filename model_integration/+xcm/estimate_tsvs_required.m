function [nt_max nt_tot nt_to nt_through Tacmat] = estimate_tsvs_required(Ng,S,k,p,alpha)
% estimate number of TSVs required in a design
% === INPUTS ===
% Ng - number of gates in the entire system
% S - number of layers in the system
% k - Rent constant of the system
% p - Rent exponent of the syystem
% alpha - input terminal fraction of the system
% === OUTPUTS ===
% nt_max - most TSVs required in a single layer
% nt_tot - total tsvs required in each layer
% nt_to - TSVs ending on each layer
% nt_through - TSVs passing through each layer

% Ng = 1e9;
% S = 8;

% k = 3.75;
% p = 0.6;
% alpha = 0.8;

% Round things out a bit so we have the same number of gates per layer
Ns = round(Ng/S);
Ng = Ns*S;


if S == 1
    nint = 0;
else
    nint = 0:S-2; % possible numbers of interstitial layers
end

Tac = 2*k*(Ns*(1+nint)).^p - k*( Ns*(2+nint)).^p - k*(nint*Ns).^p;

%% set up matrices
Tacmat = zeros(S,S); % holds number of connections between tiers
nmat = zeros(S,S);  % holds number of interstitial layers between tiers
for i=1:S-2
    for j = i+2:S
        nmat(i,j) = j - (i+1);
        nmat(j,i) = j - (i+1);
    end
end

for i=1:S
    for j=i:S
        Tacmat(i,j) = Tac(nmat(i,j)+1);
        Tacmat(j,i) = Tacmat(i,j);
    end
    Tacmat(i,i) = 0;
end

%% Calculate number of TSVs required per layer
nt_to = sum(Tacmat);
nt_tot = zeros(1,S);
for i = 1:S
    nt_tot(i) = sum(sum(Tacmat(1:i,i:end)));
end

% nt_through2 = zeros(1,S);
% for i = 2:S-1
%     p = i-1;
%     q = S-i-1;
%     nt_through2(i) = sum(sum(Tacmat(1:p,S-q:S)));
% end

nt_through = nt_tot - nt_to;

nt_req = nt_to + nt_through;
nt_max = max(nt_req);
        
    