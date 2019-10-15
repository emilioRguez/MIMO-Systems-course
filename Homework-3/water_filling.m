% Emilio Rodriguez
% Homework 3

% Waterfilling algorithm

clc; clear; close all;

% MIMO 2x2 

nt = 2;
nr = 2;

Etx = 1;

SNR = -5:15;                % [dB]
snr = 10.^(SNR/10);

Cx = (Etx/nt) * eye(nt);    % Power contraint

C = zeros(1,length(snr));

sigma_w = sqrt(Etx./snr);

% Rayleigh Channel
H = sqrt(1/2)*(randn(nr,nt) + 1i*randn(nr,nt));

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2*eye(nt);
    
    phi = real(eig(H'*inv(Cw)*H));
    step = 0.01;
    lamb = 0;
    psi = max(1/(log(2)*lamb)-1./phi,0);

    while sum(psi)-Etx > 1e-6
        lamb = lamb + step;
        psi = max(1/(log(2)*lamb)-1./phi,0);
    end
    
        C(k) = sum(log2(1 + phi.*psi));
end

% Ploting I(x,y) vs. SNR
plot(SNR,C);
title ('Channel Capacity');
xlabel('SNR');
ylabel('Capacity');
