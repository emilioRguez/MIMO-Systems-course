% Emilio Rodriguez
% Homework 6

% Implementation of ML, MMSE, ZF and MF detectors with QPSK modulation.

clear; clc; close all;

% MIMO 2x2

nt = 2;

nr = 8;

Etx = 1;

SNR = -5:5;                % [dB]
snr = 10.^(SNR/10);

sigma_w = sqrt(Etx./snr);

iter = 80000;            % Iterations number

% Alphabet with QPSK symbols
A = 1/sqrt(2) * [1+1i 1-1i -1+1i -1-1i];

% Grey codification
A_bit = [1 1; 1 0; 0 1; 0 0];% Outputs to MMSE Detector

pairs_tx = [npermutek(A,2); A(1) A(1); A(2) A(2); A(3) A(3); A(4) A(4)];

% Outputs to ML Detector
X_ML = zeros(size(pairs_tx,1),1);
BER_ML = zeros(iter,1);
BER_F_ML = zeros(length(snr),length(nr));

% Outputs to MMSE Detector
X_MMSE = zeros(length(A),2);
BER_MMSE = zeros(iter,1);
BER_F_MMSE = zeros(length(snr),length(nr));

% Outputs to ZF Detector
X_ZF = zeros(length(A),1);
BER_ZF = zeros(iter,1);
BER_F_ZF = zeros(length(snr),length(nr));

% Outputs to MF Detector
X_MF = zeros(length(A),1);
BER_MF = zeros(iter,1);
BER_F_MF = zeros(length(snr),length(nr));
           
for k = 1:length(snr)           
        
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

        %---> ML Detector
        for j = 1:size(pairs_tx,1)

            x_d = pairs_tx(j,:);

            % ML detection
            X_ML(j) = (norm(y - H * x_d.'))^2;       

        end
        
        % Dectected symbol
        [min_dist,pos] = min(X_ML);

        symb_pos_det_1 = find(abs(A-pairs_tx(pos,1)) <= 1e-10);
        symb_pos_det_2 = find(abs(A-pairs_tx(pos,2)) <= 1e-10);

        x_det_bit = [A_bit(symb_pos_det_1,:); A_bit(symb_pos_det_2,:)];

        BER_ML(l) = biterr(x_det_bit, x_bit);
        
        
        %---> MMSE Detector
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
        
        
        %---> ZF Detector
        G_ZF = (H' * H)^-1 * H';
        
        y_zf = G_ZF * y;    Cw = (sigma_w(k))^2*eye(nt);
        
        for j = 1:length(A)
            for p = 1:nt
                X_ZF(j,p) = (norm(y_zf(p) - A(j)))^2;
            end
        end
                
        % Dectected symbol
        min_dist_zf = zeros(nt,1);
        pos_zf = zeros(nt,1);
        
        for j = 1:nt
            [min_dist_zf(j),pos_zf(j)] = min(X_ZF(:,j));
        end
        
        x_det_bit_zf = [A_bit(pos_zf(1),:); A_bit(pos_zf(2),:)];

        BER_ZF(l) = biterr(x_det_bit_zf, x_bit);
        
        
        %---> MF Detector
        a = 1;
        G_MF = a * H';
        
        y_mf = G_MF * y;
        
        for j = 1:length(A)
            for p = 1:nt
                X_MF(j,p) = (norm(y_mf(p) - A(j)))^2;
            end
        end
                
        % Dectected symbol
        min_dist_mf = zeros(nt,1);
        pos_mf = zeros(nt,1);
        
        for j = 1:nt
            [min_dist_mf(j),pos_mf(j)] = min(X_MF(:,j));
        end
        
        x_det_bit_mf = [A_bit(pos_mf(1),:); A_bit(pos_mf(2),:)];

        BER_MF(l) = biterr(x_det_bit_mf, x_bit);

    end

    BER_F_ML(k) = sum(BER_ML)/(iter*4);
    BER_F_MMSE(k) = sum(BER_MMSE)/(iter*4);
    BER_F_ZF(k) = sum(BER_ZF)/(iter*4);
    BER_F_MF(k) = sum(BER_MF)/(iter*4);    

end


% Ploting
semilogy(SNR,BER_F_ML,'r-');
hold on
semilogy(SNR,BER_F_MMSE,'b-');
semilogy(SNR,BER_F_ZF,'g-');
semilogy(SNR,BER_F_MF,'m-');
legend('ML','MMSE','ZF','MF');
title ('Detectors');
xlabel('SNR');
ylabel('BER');