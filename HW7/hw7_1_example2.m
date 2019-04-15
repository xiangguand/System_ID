%{
Homework7 question 3
Symbol:
    R : correlation
    S : density spectrum
    Rnn :   Auto-correlation, noise and noise
    Rss :   Auto-correlation, sine and sine
    Rxy :   Cross-correlation, input and output
    Rs1s2 : Cross-correlation, sine_freq1 and sine_freq2, two wave by
            differenct frequency.
%}
clc, clear, close all

% Define parameters
Fs = 200;   % sampling frequency
Ts = 1/Fs;
delta_t = Ts;
L = 200;
x_t = linspace(0, Ts*L, L);

% simulation signal, two sin conbine together
sin_freq1 = 30;
sin_Am1 = 3;
sin_freq2 = 70;
sin_Am2 = 2;
t_dot = 1:1:L;

noise_mu = 0;
noise_sigma = 1;
noise = normrnd(noise_mu, noise_sigma, [1, L]);
x = sin_Am1*sin(2*pi*sin_freq1*t_dot*Ts) + sin_Am2*sin(2*pi*sin_freq2*t_dot*Ts);
y = x + noise;
figure();
plot(x_t, y);
grid on;
xlabel("Time (s)");
title("output signal");

figure();
plot(x_t, noise);
grid on;
xlabel("Time (s)");
title("Noise");

fft_y = abs(fft(y));
x_freq_raw = linspace(0, Fs/2, L/2);
figure();
plot(x_freq_raw, fft_y(1, 1:L/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Raw output frequency spectrum");

% Do circular correlation algorithm
N = int32(L/2);
R_xy = zeros([1, N]);
for i=1:N
    temp = 0;
    for j=1:N
        temp = temp + x(1,j+i)*y(1,j);
    end
    R_xy(1, i) = temp;
end
R_xy = R_xy/double(N);
S_xy = abs(fft(R_xy));

x_freq = linspace(0, Fs/2, double(N)/2);
x_time = linspace(0, Ts*double(N), double(N));
figure();
plot(x_time, R_xy);
grid on;
xlabel("Time (s)");
title("Circular correlation R_{xy}");

figure();
plot(x_freq, S_xy(1, 1:double(N)/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{xy}");


