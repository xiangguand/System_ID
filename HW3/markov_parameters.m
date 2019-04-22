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
   
    y = Y * V    % h is markov parameters
    Y: qxm
    V: mxL
%}
% Define dimension
Fs = 1;
Ts = 1/Fs;
delta_t = Ts;

L = 100;    % Create 100 terms of markov parameters

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
u0 = para_struct.input_vect;
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t);
figure();
% impulse response
impulse(sysd);

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

e = eig(Ad)

% Calculate system Markov parameters
% using symbol hn(k)
Y = zeros([q, m*L]);
U = zeros([m*L, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = Dd;
   else
       Y(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
   end
   U(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end

yn = Y*U;

t_dot = linspace(1, L*Ts, L);
figure();
plot(t_dot, yn(1, :));
xlabel('time (s)');
ylabel('output');
title('y1 Markov parameters');
grid on;

figure();
plot(t_dot, yn(2, :));
xlabel('time (s)');
ylabel('output');
title('y2 Markov parameters');
grid on;

