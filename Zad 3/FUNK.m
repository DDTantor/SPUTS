function val = FUNK(t, v, A, f, epsilon, c, i0, e1, h, n)
    val = A * v + f(v) / epsilon + [c * ones(n, 1) / epsilon; c * ones(n, 1)]...
          + 2 * epsilon / h * i0(t) * e1;
end