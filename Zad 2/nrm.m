function nr = nrm(sys_H)
    P = lyap(sys_H.A, sys_H.B * sys_H.B');
    nr = sqrt(trace(sys_H.C * P * sys_H.C'));
end