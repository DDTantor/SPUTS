function [A_r, B_r, C_r] = IRKA(A, B, C, E, red, greska, pauza, eps)
    prevD = zeros(red, 1);
    if greska
        sys_G = ss(A, B, C, 0);
        g = figure;
        movegui(g, [400, 200]);
        nrmG = nrm(sys_G);
        f = figure;
        movegui(f, [1000, 200]);
    end
    
    [n, m] = size(B);
    [p, ~] = size(C);
    
    sigma = 1e-1 + (1e3 - 1e-1) * rand(red, 1);
    r = 1e-1 + (1e3 - 1e-1) * rand(m, red);
    e = 1e-1 + (1e3 - 1e-1) * rand(p, red);
    
    V = zeros(n, red);
    W = zeros(n, red);
    lim = 100;
    residual = 1e9;
    for k = 1:lim
        fprintf('RED = %d KORAK = %d ', red, k);
        for i = 1:red
            V(:, i) = (sigma(i) * E - A) \ (B * r(:, i));
            W(:, i) = ((sigma(i) * E - A)') \ (C' * e(:, i));
        end
        
        [V, ~] = qr(V, 0);
        [W, ~] = qr(W, 0);
        
        E_p = W' * E * V;
        A_p = W' * A * V;
        B_p = W' * B;
        C_p = C * V;
        
        [X, D] = eig(A_p, E_p);
        D = diag(D);
        
        r = C_p * X;
        e = ((E_p * X) \ B_p)';
        sigma = -D;
        
        A_r = E_p \ A_p;
        B_r = E_p \ B_p;
        C_r = C_p;
        
        if greska
            figure(g);
            sys_H = ss(E_p \ A_p, E_p \ B_p, C_p, sys_G.D);
            fprintf('Relativna greska u %d-tom koraku: %f\n', k, nrm(sys_H - sys_G) / nrmG);
            bode(sys_H);
            hold on
            bode(sys_G);
            hold off
            figure(f);
            plot(real(D), imag(D), 'b+');
            hold on;
            plot(real(prevD), imag(prevD), 'r+');
            hold off;
            pause(pauza);
        end
        residual = rac_greska(prevD, D);
        fprintf('SIGMA = %e\n', residual);
        prevD = D;
        if residual < eps
            break
        end
        pause(pauza)
    end
end