function val = FUNKPODDEIM(t, v, A, g, r, II, VUPTU, V, epsilon, c, i0, h, m)
    rr = zeros(r, 1);
    for i = 1:m
        if II(i) > r
          r(i) = g(V(II(i), :) * v) / epsilon;
        else
          r(i) = c; 
        end
    end
    
    rr(1) = rr(1) + 2 * epsilon / h * i0(t);
    val = A * v + VUPTU * rr;
end