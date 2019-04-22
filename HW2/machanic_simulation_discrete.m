clc, clear, close all

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('example_machanic_n2.json');
% samping
delta_t = 0.01;

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(para_struct.state_sz))*Bc;
Cd = Cc;
Dd = Dc;

Ad
Bd
Cd
Dd
para_struct

