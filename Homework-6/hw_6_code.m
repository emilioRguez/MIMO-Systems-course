% Emilio Rodriguez
% Homework 6

% Implementation of a ML Detector with QPSK modulation.

clear; clc; close all;

% MIMO 2x2

nt = 2;

nr = [2,4,8];

Etx = 1;

SNR = -5:5:20;                % [dB]
snr = 10.^(SNR/10);

sigma_w = sqrt(Etx./snr);

iter = 10000;            % Iterations number

% Alphabet with QPSK symbols
A = 1/sqrt(2) * [1+1i 1-1i -1+1i -1-1i];

% Grey codification
A_bit = [1 1; 1 0; 0 1; 0 0];

pairs_tx = [npermutek(A,2); A(1) A(1); A(2) A(2); A(3) A(3); A(4) A(4)];

% Outputs
X_ML = zeros(size(pairs_tx,1),1);
BER = zeros(iter,1);
BER_F = zeros(length(snr),length(nr));

for p = 1:length(nr)
    
    % Rayleigh Channel
    H = sqrt(1/2)*(randn(nr(p),nt) + 1i*randn(nr(p),nt));
    
    for k = 1:length(snr)           

        for l = 1:iter

            % Gaussian Noise vector
            w = randn(nr(p),1) * sigma_w(k);

            % Transmitted signal
            pos_1 = randi(4);
            pos_2 = randi(4);

            x = [A(pos_1); A(pos_2)];

            x_bit = [A_bit(pos_1,:); A_bit(pos_2,:)];

            % Received signal
            y = H * x + w;

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

            BER(l) = biterr(x_det_bit, x_bit);

        end

        BER_F(k,p) = sum(BER)/(iter*4);

    end
end


% Ploting
semilogy(SNR,BER_F(:,1),'r-');
hold on
semilogy(SNR,BER_F(:,2),'b-');
semilogy(SNR,BER_F(:,3),'g-');
legend('nr = 2','nr = 4','nr = 8');
title ('ML Detection');
xlabel('SNR');
ylabel('BER');