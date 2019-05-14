clc, clear, close all
[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')

% analog
ssa = ss(Ac, Bc, Cc, Dc);
figure();
bode(ssa);
grid on;



