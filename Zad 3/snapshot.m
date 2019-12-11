function [A, V, W, Fs] = snapshot(n, nt, L, T, plt)
    %n = 100;
    %nt = 100;
    epsilon = 0.015;
    b = 0.5;
    gamma = 2;
    c = 0.05;
    %L = 1;
    %T = 8;
    h = L / (n - 1);
    g = @(t) t * (t-0.1) * (1-t);
    f = @(v) [arrayfun(g, v(1:n)); zeros(n, 1)];
    i0 = @(t) 50000 * t^3 * exp(-15*t);
    e1 = zeros(2 * n, 1);
    e1(1) = 1;
    dx = linspace(0, L, n);
    dt = linspace(0, T, nt);

    A11 = epsilon * 1/h^2 * (diag(ones(n-1, 1), -1)...
                           + diag(ones(n-1, 1), 1)... 
                           - 2 * eye(n));

    A11(1, 2) = epsilon * 1/h^2 * 2;
    A11(n - 1, n) = epsilon * 1/h^2 * 2;
    A12 = -eye(n) / epsilon;
    A21 = b * eye(n);
    A22 = -gamma * eye(n);
    A = [A11 A12
         A21 A22];

    F = @(t, v) FUNK(t, v, A, f, epsilon, c, i0, e1, h, n);
    
    [~, y] = ode23(F, dt, zeros(2 * n, 1));
    y = y';
    V = y(1:n, 1:nt);
    W = y(n+1:2*n, 1:nt);
    if plt == true
      f = figure();
      for i = 1:n
          plot3(dt, V(i, :), W(i, :), 'b');
          hold on;
      end
      hold off
      title('Simulacija punog sustava');
    end
    
    Fs = zeros(2 * n, nt);
    for i = 1:nt
       Fs(:, i) = [arrayfun(g, y(1:n, i)); zeros(n, 1)];
    end
end