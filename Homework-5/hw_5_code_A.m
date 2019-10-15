% Emilio Rodriguez
% Homework 5
% Item A

% Single-user MIMO
% This code needs the previous installation of the CVX toolbox. Available 
% in http://cvxr.com/cvx/download/
% This code shows that the equivalent expressions to low and high values 
% of SNR can describe the bahavior of the capacity.

clc; clear; close all;

% MIMO 2x2

nt = 2;
nr = 2;

Etx = 1;

SNR = -10:10;                % [dB]
snr = 10.^(SNR/10);

Cx = (Etx/nt) * eye(nt);    % Power contraint

C = zeros(1,length(snr));
C_high_SNR = zeros(1,length(snr));
C_low_SNR = zeros(1,length(snr));

sigma_w = sqrt(Etx./snr);

% Rayleigh fading Channel
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
    
    % Original Channel Capacity
    C(k) = log2(real(det(eye(nt) + inv(Cw)*H*Cx*H')));
    
    % Capacity equivalent to the high SNR values
    C_high_SNR(k) = min(nr,nt)*log2(Etx/sigma_w(k)^2);
    
    % Capacity equivalent to the low SNR values
    C_low_SNR(k) = nr*(Etx/sigma_w(k)^2)*log2(exp(1));
 
end


% Ploting
plot(SNR,C);
hold on
plot(SNR,C_high_SNR,'r');
hold on
plot(SNR,C_low_SNR,'g');
% plot(SNR,I);
title ('Channel Capacity');
xlabel('SNR');
ylabel('Capacity');
legend('Original Capacity','Approximation for high SNR'...
    ,'Approximation for low SNR');
