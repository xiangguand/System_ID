clc, clear, close all
%{
  @System ID Homework4
  @Please design an observer using the discrete version of A,B,C,D
  homework3 such that its G makes all the eigenvalues of /A = A + G*C are 
  zero and use the program of homework 3 to compute the observer Markov parameters.
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 100
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    A_bar matrix: nxn
    B matrix: nxm
    B_bar matrix: nx(m+q)
    G matrix: nxq
    C matrix: qxn
    D matrix: qxm
   
    y = Y * V    % Y is markov parameters
    y: qxL
    Y: qx[(q+m)(L-1)+m]
    V: [(q+m)(L-1)+m]xL
    v: (m+q)x1
%}

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
% Define dimension
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;
Fs = 10;
Ts = 1/Fs;
delta_t = Ts;
L = 500;    % Create 500 terms of markov parameters
u0 = [1;1];
y0 = [0;0];
v0 = [u0;y0];

% Design observer feedback gain
Gd = [-0.9897511 ,  0.35442173; 0.19670007, -1.2367382;-0.00490059,-0.9946962; -0.09799996,  0.09871061]

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t)

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

Ad_bar = Ad + Gd*Cd
Bd_bar = [Bd+Gd*Dd, -Gd]
eig_Ad_bar = eig(Ad_bar)

% Calculate system Markov parameters
% using symbol Yn(k)
% Discrete
Y = zeros([q, (q+m)*(L-1)+m]);
V = zeros([(q+m)*(L-1)+m, L]);
step = 1;
i = 1;
while step < (q+m)*(L-1)+m
    if step == 1
       Y(:, step:step+1) = Dd;
       V(step:step+1, i) = u0;
       step = step + q;
    else
        Y(:, step:step+3) = Cd * Ad_bar^(i-2) * Bd_bar;
        V(step:step+3, i) = v0;
        step = step + n;
    end
     i = i + 1;
end
if L ~= i-1
   display('Data length error ...');
   pause(5);
end

y = Y*V;

t_dot = linspace(1, L*Ts, L);
figure();
plot(t_dot, y(1, :));
xlabel('time (s)');
ylabel('output');
title('y1 Markov parameters (use observer)');
grid on;

figure();
plot(t_dot, y(2, :));
xlabel('time (s)');
ylabel('output');
title('y2 Markov parameters (use observer)');
grid on;

