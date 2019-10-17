% Emilio Rodriguez
% Homework 6

% Implementation of a ML Detector.

clear; clc; close all;

% MIMO 2x2

nr = 2;
nt = 2;

I_nt = eye(nt);

Etx = 1;

SNR = -5:15;                % [dB]
snr = 10.^(SNR/10);

Cx = (Etx/nt) * I_nt;    % Power contraint

C = zeros(1,length(snr));

sigma_w = sqrt(Etx./snr);


% Rayleigh Channel
H = sqrt(1/2)*(randn(nr,nt) + 1i*randn(nr,nt));

% Alphabet with QPSK symbols
A = 1/sqrt(2) * [1+1i 1-1i -1+1i -1-1i];

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2 * I_nt;
    
    % I need to compute the BER of QPSK, which depends of SNR
    
    cvx_begin sdp quiet
        variable x hermitian
        minimize()
        subject to
            Cx >= 0;
            trace(Cx) <= Etx;
    cvx_end
    
    C(k) = log2(real(det(eye(nt) + inv(Cw)*H*Cx*H')));        
 
end


% Ploting I(x,y) vs. SNR
plot(SNR,C);
title ('Channel Capacity');
xlabel('SNR');
ylabel('Capacity');