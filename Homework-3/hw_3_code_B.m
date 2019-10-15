% Emilio Rodriguez
% Homework 3
% Item B

% In this script is plotted the Chennel Capacity.

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
    
    cvx_begin sdp quiet
        variable Cx(nt,nr) hermitian
        maximize(log_det(eye(nt) + inv(Cw)*H*Cx*H'))
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
