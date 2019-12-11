function [U, S] = POD(V, r)
    [U, S, ~] = svd(V, 'econ');
    U = U(:, 1:r);
end