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
   
    y = h * V    % h is markov parameters
    y: qxL
    h: qx[(q+m)(L-1)+m]
    V: [(q+m)(L-1)+m]xL
    v: (m+q)x1
%}
[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
% Define dimension
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;
L = 100;    % Create 100 terms of markov parameters
u0 = [1;1];
y0 = [0;0];
v0 = [u0;y0];

% Design observer feedback gain
Gc = [0 0;0 0;0 0;0 0];
Gc(:,1) = -Ac(:,1);
Gc(:,2) = -Ac(:,3);

Ac_bar = Ac + Gc*Cc
Bc_bar = [Bc + Gc*Dc -Gc]
e = eig(Ac_bar)
% Calculate system Markov parameters
% using symbol hn(k)
% Continuous
h = zeros([q, (q+m)*(L-1)+m]);
V = zeros([(q+m)*(L-1)+m, L]);
yn = zeros([q, L]);
step = 1;
i = 1;
while step < (q+m)*(L-1)+m
    if step == 1
       h(:, step:step+1) = Dc;
       V(step:step+1, i) = u0;
       step = step + q;
    else
        h(:, step:step+3) = Cc * Ac_bar^(i-2) * Bc_bar;
        V(step:step+3, i) = v0;
        step = step + n;
    end
     i = i + 1;
end
if L ~= i-1
   display('Data length error ...');
   pause(5);
end

yn = h*V;

x_dot = linspace(1, L, L);
figure();
plot(x_dot, yn(1, :));
xlabel('k samples');
ylabel('output');
title('y1 Markov parameters');
grid on;

figure();
plot(x_dot, yn(2, :));
xlabel('k samples');
ylabel('output');
title('y2 Markov parameters');
grid on;

