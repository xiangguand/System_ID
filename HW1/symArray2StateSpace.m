function [T, A, B, C, D] = symArray2StateSpace(sym_Gs)
    [num, den] = numden(sym_Gs);
    num = sym2poly(num(1));
    den = sym2poly(den(1));
    T = tf(num, den);
    [A, B, C, D] = tf2ss(num, den);   % state space
end

