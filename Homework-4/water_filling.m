% Emilio Rodriguez
% Homework 4

% Waterfilling algorithm applied to Keyhole Channel.

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

% Keyhole Channel
a = randi(1,nt) + 1i*randi(1,nt);
b = randi(1,nr) + 1i*randi(1,nr);
H = (b.*a');
rank(H)

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2*eye(nt);
    
    step = 0.01;
    lamb = 0;
    
    phi = real(eig(H'*inv(Cw)*H));
   
    if lamb > phi/log(2)
        psi = 0;
    else
        psi = max(1/(log(2)*lamb)-1./phi,0);
    end
    

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
