% Emilio Rodriguez
% Homework 8

% Implementation of standard MMSE, MMSE 1 bit and improved MMSE 1 bit
% detectors with QPSK modulation.

clear; clc; close all;

% MIMO 2x2

nt = 2;

nr = 8;

Etx = 1;

SNR = -5:5;                % [dB]
snr = 10.^(SNR/10);

sigma_w = sqrt(Etx./snr);

Cx = Etx/nt * eye(nt);

iter = 800;            % Iterations number

% Alphabet with QPSK symbols
A = 1/sqrt(2) * [1+1i 1-1i -1+1i -1-1i];

% Grey codification
A_bit = [1 1; 1 0; 0 1; 0 0]; % Outputs to MMSE Detector

% pairs_tx = [npermutek(A,2); A(1) A(1); A(2) A(2); A(3) A(3); A(4) A(4)];

% Outputs to MMSE Detector
X_MMSE = zeros(length(A),1);
BER_MMSE = zeros(iter,1);
BER_F_MMSE = zeros(length(snr),length(nr));

% Outputs to MMSE Detector with 1 bit ADCs
X_MMSE_1_bit = zeros(length(A),1);
BER_MMSE_1_bit = zeros(iter,1);
BER_F_MMSE_1_bit = zeros(length(snr),length(nr));

% Outputs to Normal MMSE Detector with 1 bit ADCs
X_MMSE_1_bit_norm = zeros(length(A),1);
BER_MMSE_1_bit_norm = zeros(iter,1);
BER_F_MMSE_1_bit_norm = zeros(length(snr),length(nr));
           
for k = 1:length(snr)           
        
    Cw = (sigma_w(k))^2 * eye(nr);
    
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

        % Received signal
        y = H * x + w;
        
        %---> Standard MMSE Detector
        G_MMSE = (nt * (sigma_w(k))^2 * eye(nt)/Etx + H' * H)^-1 * H';
        
        y_mmse = G_MMSE * y;
        
        for j = 1:length(A)
            for p = 1:nt
                X_MMSE(j,p) = (norm(y_mmse(p) - A(j)))^2;
            end
        end
                
        % Dectected symbol
        min_dist_mmse = zeros(nt,1);
        pos_mmse = zeros(nt,1);
        
        for j = 1:nt
            [min_dist_mmse(j),pos_mmse(j)] = min(X_MMSE(:,j));
        end
        
        x_det_bit_mmse = [A_bit(pos_mmse(1),:); A_bit(pos_mmse(2),:)];

        BER_MMSE(l) = biterr(x_det_bit_mmse, x_bit);
        
        
        %---> MMSE Detector with 1-bit ADCs
        
         G_MMSE_1_bit = sqrt(pi/2)*Cx*H' * (diag(diag(H*Cx*H'+Cw)))^(-1/2)...
             * (asin((diag(diag(H*Cx*H'+Cw)))^(-1/2)* real(H*Cx*H'+Cw)*...
             (diag(diag(H*Cx*H'+Cw)))^(-1/2)) + 1i*asin((diag(diag(H*...
             Cx*H'+Cw)))^(-1/2)* imag(H*Cx*H'+Cw)*(diag(diag(H*Cx*H'...
             +Cw)))^(-1/2)))^(-1);
         
         y_mmse_1_bit = G_MMSE_1_bit * y;
                
        for j = 1:length(A)
            for p = 1:nt
                X_MMSE_1_bit(j,p) = (norm(y_mmse_1_bit(p) - A(j)))^2;
            end
        end
                
        % Dectected symbol
        min_dist_mmse_1_bit = zeros(nt,1);
        pos_mmse_1_bit = zeros(nt,1);
        
        for j = 1:nt
            [min_dist_mmse_1_bit(j),pos_mmse_1_bit(j)] = min(X_MMSE_1_bit(:,j));
        end
        
        x_det_bit_mmse_1_bit = [A_bit(pos_mmse_1_bit(1),:); A_bit(pos_mmse_1_bit(2),:)];

        BER_MMSE_1_bit(l) = biterr(x_det_bit_mmse_1_bit, x_bit);
        
        
        %---> MMSE Detector with 1-bit ADCs
             
        y_1bit = sign(real(H*x + w)) + 1i*sign(imag(H*x + w));
        
        y_mmse_1bit_norm = G_MMSE * y_1bit;
        
        for j = 1:length(A)
            for p = 1:nt
                X_MMSE_1_bit_norm(j,p) = (norm(y_mmse_1bit_norm(p) - A(j)))^2;
            end
        end
                
        % Dectected symbol
        min_dist_mmse_1_bit_norm = zeros(nt,1);
        pos_mmse_1_bit_norm = zeros(nt,1);
        
        for j = 1:nt
            [min_dist_mmse_1_bit_norm(j),pos_mmse_1_bit_norm(j)] = min(X_MMSE_1_bit_norm(:,j));
        end
        
        x_det_bit_mmse_1_bit_norm = [A_bit(pos_mmse_1_bit_norm(1),:); A_bit(pos_mmse_1_bit_norm(2),:)];

        BER_MMSE_1_bit_norm(l) = biterr(x_det_bit_mmse_1_bit_norm, x_bit);
        
    end

    BER_F_MMSE(k) = sum(BER_MMSE)/(iter*4);
    BER_F_MMSE_1_bit(k) = sum(BER_MMSE_1_bit)/(iter*4);
    BER_F_MMSE_1_bit_norm(k) = sum(BER_MMSE_1_bit_norm)/(iter*4);

end


% Ploting
semilogy(SNR,BER_F_MMSE,'b-');
hold on
semilogy(SNR,BER_F_MMSE_1_bit,'g-');
semilogy(SNR,BER_F_MMSE_1_bit_norm,'m-');
legend('Standard MMSE','MMSE with 1-bit','Normal MMSE with 1-bit');
title ('Detectors');
xlabel('SNR');
ylabel('BER');