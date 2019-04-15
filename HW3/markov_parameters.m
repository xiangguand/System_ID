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
delta_t = 0.01;

L = 2000;    % Create 100 terms of markov parameters

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
u0 = para_struct.input_vect;
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t);

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

e = eig(Ad)

% Calculate system Markov parameters
% using symbol hn(k)
h = zeros([q, m*L]);
U = zeros([m*L, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       h(:, step:step+m-1) = Dd;
   else
       h(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
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

