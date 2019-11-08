% Emilio Rodriguez
% Homework 9

% Implementation of the linear precoding with MMSE criterion using
% QPSK modulation.

clear; clc; close all;

% MIMO 4x2

nt = 4;

nr = 2;

Etx = 1;

SNR = -5:15;                 % [dB]
snr = 10.^(SNR/10);

sigma_w = sqrt(Etx./snr);

Int = eye(nt);

Inr = eye(nr);

Cx = Inr;

iter = 8000;                 % Iterations number

% Alphabet with QPSK symbols
A = 1/sqrt(2) * [1+1i 1-1i -1+1i -1-1i];

% Grey codification
A_bit = [1 1; 1 0; 0 1; 0 0];

% Outputs to Detector
X_det = zeros(length(A),1);
BER = zeros(iter,1);
BER_F = zeros(length(snr),length(nr));

    
for k = 1:length(snr)           

    Cw = (sigma_w(k))^2 * Inr;

    for l = 1:iter

        % Rayleigh Channel
        H = sqrt(1/2)*(randn(nr,nt) + 1i*randn(nr,nt));

        % Gaussian Noise vector
        w = (randn(nr,1) + 1i*randn(nr,1)) * sqrt(1/2)*sigma_w(k);

        % Transmitted signal
        pos_1 = randi(4);
        pos_2 = randi(4);

        x = [A(pos_1); A(pos_2)];

        x_bit = [A_bit(pos_1,:); A_bit(pos_2,:)];

        % Parameters of Precoding
        G = Inr;
        lamb_f2 = trace(G*Cw*G')/Etx;
        f = sqrt(trace((H'*G'*G*H + trace(G*Cw*G')/Etx*Int)^-2 * H'*G'...
            *Cx*G*H)/Etx);
                
        % Precoding
        P_MMSE = 1/f * (H'*G'*G*H + lamb_f2 * Int)^-1 * H'*G';
        
        % Received signal
        y = H * P_MMSE * x + w;

        %---> Detector
        y_det = G * y;

        for j = 1:length(A)
            for p = 1:nr
                X_det(j,p) = (norm(y_det(p) - A(j)))^2;
            end
        end

        % Dectected symbol
        min_dist = zeros(nr,1);
        pos = zeros(nr,1);

        for j = 1:nr
            [min_dist(j),pos(j)] = min(X_det(:,j));
        end

        x_det_bit = [A_bit(pos(1),:); A_bit(pos(2),:)];

        BER(l) = biterr(x_det_bit, x_bit);

    end

    BER_F(k) = sum(BER)/(iter*4);

end


% Ploting
semilogy(SNR,BER_F,'LineWidth',1.3);
title ('Detectors with precoding');
xlabel('SNR');
ylabel('BER');