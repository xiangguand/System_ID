%{
This file is HW7 question 1
Do homework7 question1, to calculate Circular Correlation, R_{xy}, of two
sequences, x[k], y[k]
xcorr
%}
clc, clear, close all
% define parameters
Fs = 300;
Ts = 1/Fs;
delta_t = Ts;
L = 100;
t_dot = 1:1:L;
noise1 = rand([1, L])*3;
noise2 = rand([1, L])*3;
sin_freq1 = 15;
sin_freq2 = 30;
A1 = 5;
A2 = 3;
x = A1*sin(2*pi*sin_freq1*t_dot*delta_t) + noise1;
y = A2*sin(2*pi*sin_freq2*t_dot*delta_t) + noise2;

figure();
x_t = linspace(0, L/Fs, L);
hold on;
grid on;
plot(x_t, x, "r");
plot(x_t, y, "b");
legend("input x", "output y");
title("raw input and output");
xlabel("time (s)");
ylabel("amplitude");

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

x2_t = linspace(0, double(N)/Fs, N);
figure();
grid on; hold on;
plot(x2_t, R_xy(1, :));
plot(x2_t, x(1, 1:N));
plot(x2_t, y(1, 1:N));
legend("correlation R_{xy}", "input x", "output y");
title("correlation R_{xy} and raw input");
xlabel("time (s)");
ylabel("amplitude");

% note
% As we can see. After we do correlation, it can press the noise and let
% the signal be prettier.
% Warning !!! phase will change

