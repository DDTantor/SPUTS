function simulacijaPODDEIM(V, II, A, UPTU, r, m, nt, n, T, L, plt)
  epsilon = 0.015;
  b = 0.5;
  gamma = 2;
  c = 0.05;
  
  A = V' * A * V;
  g = @(t) t * (t-0.1) * (1-t);
  
  #f = @(v) [arrayfun(g, v(1:n)); zeros(n, 1)];
  i0 = @(t) 50000 * t^3 * exp(-15*t);
  dt = linspace(0, T, nt);
  h = L / (n - 1);
  VUPTU = V' * UPTU;
  F = @(t, v) FUNKPODDEIM(t, v, A, g, r, II, VUPTU, V, epsilon, c, i0, h, m);
  
  [~, y] = ode23(F, dt, zeros(2 * r, 1));
  y = y';
  V = y(1:r, :);
  W = y(r+1:2*r, :);
  if plt == true
      f = figure();
      for i = 1:r
          plot3(dt, V(i, :), W(i, :), 'b');
          hold on;
      end
      hold off
      title('Simulacija POD + DEIM');
    end
end