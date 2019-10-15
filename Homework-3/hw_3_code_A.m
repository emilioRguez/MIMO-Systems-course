% Emilio Rodriguez
% Homework 3
% Item A

% In this script is plotted the Mutual Information.

clc; clear; close all;

% MIMO 2x2

nt = 2;
nr = 2;

Etx = 1;

SNR = -5:15;            % [dB]
snr = 10.^(SNR/10);

Cx = (Etx/nt) * eye(nt);    % Power constraint

I = zeros(1,length(snr));

sigma_w = sqrt(Etx./snr);

% Rayleigh Channel

H = sqrt(1/2)*(randn(nr,nt) + 1i*randn(nr,nt));

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2*eye(nt); 
        
    I(k) = log2(real(det(eye(nt) + H'*inv(Cw)*H*Cx)));
end


% Ploting I(x,y) vs. SNR
plot(SNR,I);
title ('Mutual Information');
xlabel('SNR');
ylabel('I(x,y)');
