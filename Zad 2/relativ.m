function res = relativ(sys_G, sys_Gt)
    sys_H = sys_Gt - sys_G;
    P = lyap(sys_H.A, sys_H.B * sys_H.B');
    res = sqrt(trace(sys_H.C * P * sys_H.C'));
end