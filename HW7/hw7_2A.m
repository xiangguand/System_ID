clc, clear, close all
%{
  @System ID Homework7 question 2
  @Please write a program which can transform the representation of a
  system from state space model to the linear defference model.
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 1000
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    B matrix: nxm
    C matrix: qxn
    D matrix: qxm
    
    u: mxL
    u[k] = sin(w*k*delta_t)
    Yi: qxm
    y: qxL    

%}


[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
% Define dimension
u0 = para_struct.input_vect;
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;
L = 400;
% Decide the sampling rate according to the analog system
Fs = 2;             % Sampling frequency    
SYSTEM_DIM = 4;             
Ts = 1/Fs;             % Sampling period
delta_t = Ts;
t = (0:L-1)*Ts;        % Time vector

% Degitalize the system with the above sampling rate to obtain
% its digital model
Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;
ssd = ss(Ad, Bd, Cd, Dd);

% Calculate weighted sequence (Markov parameters)
Yi = zeros([q, m*L]);
U = zeros([m*L, L]);    % impulse input
step = 1;
i = 1;
while step < m*L
   if step == 1
       Yi(:, step:step+m-1) = Dd;
   else
       Yi(:, step:step+m-1) = Cd*(Ad^(i-2))*Bd;
   end
   U(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end

y = Yi*U;   % impulse response

Yi1 = y(1, :);
Yi2 = y(2, :);
x_time = linspace(0, L*Ts, L);   % Yi dimension is 2*m*L
figure();
hold on; grid on;
title("Weighted sequence in time domain");
plot(x_time, Yi1, "b");
plot(x_time, Yi2, "r");
xlabel("Time (s)");
ylabel("Amplitude");
legend("Yi1", "Yi2");

% HW7 question2 A
Yi1_fft = abs(fft(y(1,:)));
Yi2_fft = abs(fft(y(2,:)));
x_freq = linspace(0, Fs/2, L/2);
Yi_Am1 = Yi1_fft(1, 1:L/2);
Yi_Am2 = Yi2_fft(1, 1:L/2);
figure();
hold on; grid on;
title("Weighted sequence after FFT");
plot(x_freq, Yi_Am1, "b");
plot(x_freq, Yi_Am2, "r");
xlabel("Frequency (Hz)");
ylabel("Amplitude");
legend("Yi1_fft", "Yi2_fft");









