% Emilio Rodriguez
% Homework 5
% Item B

% Multiuser MIMO 
% Uplink configuration.
% This code needs the previous installation of the CVX toolbox. Available 
% in http://cvxr.com/cvx/download/

clc; clear; close all;

nr = 2;
nt_1 = 2;
nt_2 = 2;

K = 2;
Inr = eye(nr);
Int = eye(nt_1);

Etx_1 = 1;
Etx_2 = 1;

sigma_w = sqrt(Etx_1);


% Rayleigh fading Channel
H_1 = sqrt(1/2)*(randn(nr,nt_1) + 1i*randn(nr,nt_1));
H_2 = sqrt(1/2)*(randn(nr,nt_2) + 1i*randn(nr,nt_2));

Cw = sigma_w^2*Int;

%Point A
cvx_begin sdp quiet
    variable Cx_1(nr,nt_1) hermitian
    maximize(log_det(Inr + Cw^-1*H_1*Cx_1*H_1'))
    subject to
        Cx_1 >= 0;
        trace(Cx_1) <= Etx_1;
cvx_end

R1_A = log2(real(det(Inr + Cw^-1*H_1*Cx_1*H_1')));
    
cvx_begin sdp quiet
    variable Cx_2(nr,nt_2) hermitian
    maximize(log_det(Inr + (Cw^-1*H_2*Cx_2*H_2')*(Inr + Cw^-1*H_1*Cx_1*H_1')^-1))
    subject to
        Cx_2 >= 0;
        trace(Cx_2) <= Etx_2;
cvx_end

R2_A = log2(real(det(Inr + (Cw^-1*H_2*Cx_2*H_2')*(Inr + Cw^-1*H_1*Cx_1*H_1')^-1)));


%Point B
cvx_begin sdp quiet
    variable Cx_2(nr,nt_2) hermitian
    maximize(log_det(Inr + Cw^-1*H_2*Cx_2*H_2'))
    subject to
        Cx_2 >= 0;
        trace(Cx_2) <= Etx_2;
cvx_end

R2_B = log2(real(det(Inr + Cw^-1*H_2*Cx_2*H_2')));
    
cvx_begin sdp quiet
    variable Cx_1(nr,nt_1) hermitian
    maximize(log_det(Inr + (Cw^-1*H_1*Cx_1*H_1')*(Inr + Cw^-1*H_2*Cx_2*H_2')^-1))
    subject to
        Cx_1 >= 0;
        trace(Cx_1) <= Etx_2;
cvx_end

R1_B = log2(real(det(Inr + (Cw^-1*H_1*Cx_2*H_1')*(Inr + Cw^-1*H_2*Cx_2*H_2')^-1)));

R_A = [R1_A R2_A];
R_B = [R1_B R2_B];

title ('Channel Capacity');
plot(R_A(1),R_A(2),'.','MarkerSize',30,'Color','blue');
pl_A = line([R_A(1) R_A(1)],[R_A(2) 0],'Color','blue','LineWidth',1.3);
hold on
plot(R_B(1),R_B(2),'.','MarkerSize',30,'Color','red');
pl_B = line([0 R_B(1)],[R_B(2) R_B(2)],'Color','red','LineWidth',1.3);
xlabel('R_1');
ylabel('R_2');
