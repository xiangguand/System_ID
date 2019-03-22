clc, clear, close all
%{
  @System ID Homework3
  @Please program and compute Markov parameters.
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 100
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    B matrix: nxm
    C matrix: qxn
    D matrix: qxm
   
    y = h * V    % h is markov parameters
    h: qxm
    V: mxL
%}
% Define dimension
m = 2;
q = 2;
n = 4;
L = 100;    % Create 100 terms of markov parameters
u0 = [1;1];

% Define parameters
M = 1;   % mass
K = 1;   % spring
zeta = 0.1; % damping

Ac = [0 1 0 0;-(2*K/M) -(2*zeta/M) (K/M) (zeta/M);0 0 0 1;(K/M) (zeta/M) -(K/M) -(zeta/M)]
Bc = [0 0;1/M 0;0 0;0 1/M]
Cc = [1 0 0 0;0 0 1 0]
Dc = [0 0;0 0]
model = ss(Ac, Bc, Cc, Dc);

e = eig(Ac)

% Calculate system Markov parameters
% using symbol hn(k)
h = zeros([q, m*L]);
U = zeros([m*L, L]);
yn = zeros([q, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       h(:, step:step+m-1) = Dc;
   else
       h(:, step:step+m-1) = Cc*(Ac^(step-2))*Bc;
   end
   U(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end

yn = h*U;

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


