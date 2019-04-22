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
Fs = 10;
Ts = 1/Fs;
delta_t = Ts;

L = 2000;    % Create 2000 terms of markov parameters

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

% Calculate system Markov parameters (symbol h)
% using symbol hn(k)
h = zeros([q, m*L]);
% Random noise generates, to exite system
U = zeros([m*L, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       h(:, step:step+m-1) = Dd;
   else
       h(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
   end
   U(step:step+m-1, i) = rand([2, 1]);
   i = i + 1;
   step = step + m;
end

% output
y = h*U;
t_dot = linspace(0, L*Ts, L);

figure();
subplot(2,1,1);
plot(t_dot, y(1,:));
title("output1 random exitation");
xlabel("time (s)");
ylabel("Value");
grid on;

subplot(2,1,2);
plot(t_dot, y(2,:));
title("output2 random exitation");
xlabel("time (s)");
ylabel("Value");
grid on;


y1_fft = abs(fft(y(1, :)));
y2_fft = abs(fft(y(2, :)));
y1_fft = y1_fft(1:L/2);
y2_fft = y2_fft(1:L/2);

freq_dot = linspace(0, Fs/2, L/2);
figure();
subplot(2,1,1);
plot(freq_dot, y1_fft);
title("output1 random exitation");
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

subplot(2,1,2);
plot(freq_dot, y2_fft);
title("output2 random exitation");
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;






