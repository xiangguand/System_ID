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
Fs = 300;   % sampling frequency
Ts = 1/Fs;
delta_t = Ts;
L = 500;
x_t = linspace(0, Ts*L, L);

% simulate noise
noise_mu1 = 0;
noise_sigma1 = 1;
noise_mu2 = 0;
noise_sigma2 = 1;
noise1 = normrnd(noise_mu1, noise_sigma1, [1, L]);
noise2 = normrnd(noise_mu2, noise_sigma2, [1, L]);

% simulate signal
sin_freq1 = 10;
sin_Am1 = 10;
sin_freq2 = 100;
sin_Am2 = 5;
t_dot = 1:1:L;
sin1 = sin_Am1*sin(2*pi*sin_freq1*t_dot*Ts) + noise1;
sin2 = sin_Am2*sin(2*pi*sin_freq2*t_dot*Ts) + noise1;

figure();
plot(x_t, sin1);
grid on;
xlim([0, L*Ts]);
xlabel("Time (s)");
title("sine wave number 1");

figure();
plot(x_t, sin2);
grid on;
xlim([0, L*Ts]);
xlabel("Time (s)");
title("sine wave number 2");

figure();
plot(x_t, noise1);
grid on;
xlim([0, L*Ts]);
xlabel("Time (s)");
title("normal noise number 1");

figure();
plot(x_t, noise2);
grid on;
xlim([0, L*Ts]);
xlabel("Time (s)");
title("normal noise number 2");

% Do circular correlation algorithm
N = double(int32(L/2));
Rnn = zeros([1, N]);
Rss = zeros([1, N]);
Rn1n2 = zeros([1, N]);
Rns = zeros([1, N]);
Rs1s2 = zeros([1, N]);

for i=1:N
    temp_nn = 0;
    temp_ss = 0;
    temp_n1n2 = 0;
    temp_ns = 0;
    temp_s1s2 = 0;
    for j=1:N
        temp_nn = temp_nn + noise1(1,j+i) * noise1(1,j);
        temp_ss = temp_ss + sin1(1,j+i) * sin1(1,j);
        temp_n1n2 = temp_n1n2 + noise1(1,j+i) * noise2(1,j);
        temp_ns = temp_ns + noise1(1,j+i) * sin1(1,j);
        temp_s1s2 = temp_s1s2 + sin1(1,j+i) * sin2(1,j);
    end
    Rnn(1, i) = temp_nn / N;
    Rss(1, i) = temp_ss / N;
    Rn1n2(1, i) = temp_n1n2 / N;
    Rns(1, i) = temp_ns / N;
    Rs1s2(1, i) = temp_s1s2 / N;
end

Snn = abs(fft(Rnn));
Sss = abs(fft(Rss));
Sn1n2 = abs(fft(Rn1n2));
Sns = abs(fft(Rns));
Ss1s2 = abs(fft(Rs1s2));

x_freq = linspace(0, Fs/2, N/2);
x_time = linspace(0, Ts*N, N);

% ################################################## %
% 3.A
figure();
subplot(2,1,1);
plot(x_time, Rnn);
grid on;
xlim([0, N*Ts]);
xlabel("Time (s)");
title("Circular correlation R_{nn}");
subplot(2,1,2);
plot(x_freq, Snn(1, 1:N/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{nn}");

% ################################################## %
% 3.B
figure();
subplot(2,1,1);
plot(x_time, Rss);
grid on;
xlim([0, N*Ts]);
xlabel("Time (s)");
title("Circular correlation R_{ss}");
subplot(2,1,2);
plot(x_freq, Sss(1, 1:N/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{ss}");

% ################################################## %
% 3.C
figure();
subplot(2,1,1);
plot(x_time, Rn1n2);
grid on;
xlim([0, N*Ts]);
xlabel("Time (s)");
title("Circular correlation R_{n1n2}");
subplot(2,1,2);
plot(x_freq, Sn1n2(1, 1:N/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{n1n2}");

% ################################################## %
% 3.D
figure();
subplot(2,1,1);
plot(x_time, Rns);
grid on;
xlim([0, N*Ts]);
xlabel("Time (s)");
title("Circular correlation R_{ns}");
subplot(2,1,2);
plot(x_freq, Sns(1, 1:N/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{ns}");

% ################################################## %
% 3.E
figure();
subplot(2,1,1);
plot(x_time, Rs1s2);
grid on;
xlim([0, N*Ts]);
xlabel("Time (s)");
title("Circular correlation R_{s1s2}");
subplot(2,1,2);
plot(x_freq, Ss1s2(1, 1:N/2));
grid on;
xlabel("Frequency (Hz)");
ylabel("Amplitude");
title("Circular correlation density spectrum S_{s1s2}");


