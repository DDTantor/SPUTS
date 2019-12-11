function val = FUNKPOD(t, v, A, V, f, r, epsilon, c, i0, e1, h, n)
    cc = [c * ones(n, 1) / epsilon; c * ones(n, 1)];
    i0t = 2 * epsilon / h * i0(t) * e1;
    val = A * v + V' * (f(V * v) / epsilon + cc + i0t);
end