% Emilio Rodriguez
% Homework 2: Gemetrical Channel Model for line-of-sigh MIMO

clc; clear; close all;

nt = 4;
nr = 4;
D = 0.30; % 30 cm
f = 60e9; % 60 GHz
c = 3e8;
lamb = c/f;

mi = [0 1]; % tx
ni = [0 1];
mj = [0 1]; % rx
nj = [0 1];

dv = 0:0.0001:0.05;
dh = dv;

H = zeros(nt,nr);
P = zeros(1,length(dv));

var_x = 1/4;
var_w = 1;

for k = 1:length(dv)
    for t = 1:nt
        for r = 1:nr
            t_b = de2bi(t-1,2);
            r_b = de2bi(r-1,2);
            H(r,t) = exp(1i*((2*pi)/lamb)*(((dv(k)^2)/(2*D))*(r_b(1)-t_b(1))^2 ...
            + ((dh(k)^2)/(2*D))*(r_b(2)-t_b(2))^2));
        end
    end
    P(1,k) = real(det(eye(nt) + (var_x/var_w)*(H')*H));
end

I = log10(P)/log10(2);

% Ploting I(x,y) vs. dv
plot(dv,I);
title ('Mutual Information vs. Horizontal distance between antennas');
xlabel('dv(m)');
ylabel('I(x,y)');
