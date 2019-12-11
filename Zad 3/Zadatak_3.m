n = 50;
nt = 30;
L = 1;
T = 8;
r = 20; % POD redukcija
m = 20; % DEIM
plt = true; % true ako zelimo da se plotaju plot3

%% snapshot(diskretizacija po x, 
%%          diskretizacija po t
%%          duljina po x
%%          vrijeme po t
%%          bool - plot3 

prev = time();
[A, V, W, Fs] = snapshot(n, nt, L, T, plt);
fprintf('Vrijeme potrebno za punu simulaciju: %f\n', time() - prev);
[Uv, Sv] = POD(V, r);
[Uw, Sw] = POD(W, r);
[Uf, Sf, ~] = svd(Fs, 'econ');

if plt == true
  f = figure();
  semilogy(diag(Sv) + 1e-18, 'ro'); %% Moram dodat 1e-18 jer su neke -0
  hold on
  semilogy(diag(Sw) + 1e-18, 'g*');
  semilogy(diag(Sf) + 1e-18, 'bx');
  title('Singularne vrijednosti');
  legend('V', 'W', 'f');
  hold off
end

[II, PTU] = DEIM(Uf, m);
VV = blkdiag(Uv, Uw);
prev = time();
simulacijaPOD(VV, A, r, nt, n, T, L, plt);
fprintf('Vrijeme potrebno za POD simulaciju: %f\n', time() - prev);
prev = time();
simulacijaPODDEIM(VV, II, A, Uf(:, 1:m) * PTU, r, m, nt, n, T, L, plt);
fprintf('Vrijeme potrebno za POD + DEIM simulaciju: %f\n', time() - prev);