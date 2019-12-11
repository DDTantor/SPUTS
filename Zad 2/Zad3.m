clc
clear

A = mmread('rail_20209.mtx');
B = mmread('rail_20209_B.mtx');
C = mmread('rail_20209_C.mtx');
E = mmread('rail_20209_E.mtx');
D = zeros(1);

% IRKA(A, B, C, E, red, greska, pauza, eps)
[A_r, B_r, C_r] = IRKA(A, B(:, 6), C(2, :), E, 10, false, 0.005,1e-2);
f = figure;
sys_Gr = ss(A_r, B_r, C_r, D);
[mag, ph, w] = bode(sys_Gr, logspace(-9, 2, 200));

%%
g = figure;
magdb = exp(squeeze(mag) / 20);
semilogx(squeeze(w), 20 * log10(squeeze(mag)));