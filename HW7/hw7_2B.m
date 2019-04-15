%{
This file is HW7 question 2 B
Discuss the effect of DFT length, 10,10.25.10.5,100,100.25,100.5 period

%}
clc, clear, close all

Fs = 0.5;  % sampling frequency
Ts = 1/Fs;
L = 1000;

sin_period = [10, 10.25, 10.5, 100, 100.25, 100.5];
sin_freq = 1./sin_period
sin_dim = 6;
sin_dot = ones(sin_dim, L);
dot = 1:1:L;
A = 2;
noise = rand([1, L]);

u = zeros([sin_dim, L]);
fft_u = zeros([sin_dim, L]);
for i=1:sin_dim
   u(i, :) = A*sin(2*pi*(1/sin_period(1, i))*dot*Ts) + noise; 
   fft_u(i, :) = abs(fft(u(i, :)));
end


figure();
hold on;
grid on;
xlabel("Frequency");
ylabel("Amplitude");
xlim([0, max(sin_freq)+0.05]);
x_freq = linspace(0, Fs/2, L/2);
for i=1:sin_dim
   plot(x_freq, fft_u(i, 1:L/2));
end
legend(num2str(sin_freq(1)), num2str(sin_freq(2)), num2str(sin_freq(3)),...
       num2str(sin_freq(4)), num2str(sin_freq(5)), num2str(sin_freq(6)));


% Frequency is too small, three frequency (0.01¡B0.0099¡B0.0099), they
% overlay.

