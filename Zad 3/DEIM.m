function [II, PTU] = DEIM(U, m)
  [~, II(1)] = max(abs(U(:, 1)));
  PTU = U(II, 1);
  for l = 2:m
    c = U(II, l) \ PTU;
    r = U(:, l) - U(:, 1:(l - 1)) * c';
    [~, II(l)] = max(abs(r));
    II;
    PTU;
    U(II(l), 1:l);
    PTU = [PTU, zeros(l - 1, 1);
           U(II(l), 1:l)];
  end
end