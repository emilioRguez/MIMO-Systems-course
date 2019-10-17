% Emilio Rodriguez
% Homework 6

% Implementation of a ML Detector.

clear; clc; close all;

% MIMO 2x2

nr = 2;
nt = 2;

sigma_w = 1;

I_nr = eye(nr);

Cw = sigma_w^2*I_nr;

Cx = (Etx/nt)*I_nr;