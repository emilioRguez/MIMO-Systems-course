% Emilio Rodriguez 1913153
% Homework 3
% Item A

% In this script is plotted the Mutual Information but now was used a
% Keyhole Channel.

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

% Keyhole Channel
a = randi(1,nt) + 1i*randi(1,nt);
b = randi(1,nr) + 1i*randi(1,nr);
H = b.*a';

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2*eye(nt); 
        
    I(k) = log2(real(det(eye(nt) + H'*inv(Cw)*H*Cx)));
end


% Ploting I(x,y) vs. SNR
plot(SNR,I);
title ('Mutual Information');
xlabel('SNR');
ylabel('I(x,y)');
