clc
clear
format short e
load CDplayer.mat
D = zeros(1);
n = length(A);
B = B(:, 2);
C = C(1, :);
sys_G = ss(A, B, C, D);
u = zeros(37, 1);
for r = 3:40
    sys_Gt = balred(sys_G, r);
    u(r - 2) = nrm(sys_G - sys_Gt) / nrm(sys_G);
end

%figure;
%semilogy(u);
% A B C E red racunanje_greske(true ili false) duljina pauze
%[A_r, B_r, C_r] = IRKA(A, B, C, eye(n), 15, true, 0.005);
%sys_Gr = ss(A_r, B_r, C_r, D);
%sys_G = ss(A, B, C, D);

%fprintf('Greška: %f\n', nrm(sys_G - sys_Gr) / nrm(sys_G));

%%
for r = 3:40
    [A_r, B_r, C_r] = IRKA(A, B, C, eye(n), r, false, 0, 1e-3);
    sys_Gt = ss(A_r, B_r, C_r, D);
    u(r - 2) = nrm(sys_G - sys_Gt) / nrm(sys_G);
    fprintf('RED = %d GRESKA = %e\n', r, u(r - 2));
end

figure();
semilogy(real(u));