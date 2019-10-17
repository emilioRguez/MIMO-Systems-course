% Emilio Rodriguez
% Homework 6

% Implementation of a ML Detector with QPSK modulation.

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

pairs_tx = [npermutek(A,2); A(1) A(1); A(2) A(2); A(3) A(3); A(4) A(4)];

% Outputs
X_ML = zeros(size(pairs_tx,1),1);

for k = 1:length(snr)
    
    Cw = (sigma_w(k))^2 * I_nt;
    
    % Gaussian Noise vector
    w = randn(nr,1) * sigma_w(k);
    
    for i = 1:size(pairs_tx,1)
       
        x = pairs_tx(i);
        
        % Receiver signal
        y = H * x' + w;
        
        % ML detection
        X_ML(i) = (norm(y - H * x'))^2;
       
    end
    
    % Dectected symbol
    [x_d,pos] = min(X_ML);
    
    
    
    
    
    
   
    
    
            
 
end


% Ploting I(x,y) vs. SNR
% plot(SNR,C);
% title ('Channel Capacity');
% xlabel('SNR');
% ylabel('Capacity');