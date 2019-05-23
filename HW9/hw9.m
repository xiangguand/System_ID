clc, clear, close all
%{
  @System ID Homework9
  @Recuisive least square (RLS)
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    arn : AR model parameters numbers
    parameters A:arn x q
    gain K:arn x 1
    output y:1xq
    error e:1xq
    input x[t]:1 x arn
    P: arn x arn
    sigma: 1x1
%}
% Define dimension
Fs = 500;
Ts = 1/Fs;
m = 1;
q = 1;
AR_PARAMETERS_NUMBERS = 2;  % tune two coefficient;
delta_t = Ts;
lambda = 0.95;
step_times = 1000;
L = 2000;    % Create 2000 terms of markov parameters

% simulation
k_dots = 1:1:L;
sin_freq1 = 25;
sin_A1 = 30;
sin_freq2 = 20;
sin_A2 = 5;
noise = rand(1, L)*0.01;
y_raw = sin_A1*sin(2*pi*sin_freq1*k_dots*Ts);
SNR = 20 * log(mean((y_raw./noise).^2)^0.5)
y_raw = y_raw + noise;
figure();
plot(k_dots*Ts, y_raw);
title('simulate unknow signal')
xlabel('time (t)');
ylabel('amplitude');
grid on;

% step three
% y[t] = ar1 * y[t-1] + ar2 * y[t-2]
Y = y_raw(:, AR_PARAMETERS_NUMBERS+1:end)';
X = zeros([L-AR_PARAMETERS_NUMBERS, AR_PARAMETERS_NUMBERS]);
for i = 1:AR_PARAMETERS_NUMBERS
   X(:, i) = y_raw(:, AR_PARAMETERS_NUMBERS-i+1:L-i)';
end

% update the parameters
A = rand([AR_PARAMETERS_NUMBERS, 1]);    % initial A parameters
error = zeros([L-2*AR_PARAMETERS_NUMBERS, 1]);
i = 1;
for step = AR_PARAMETERS_NUMBERS:L-AR_PARAMETERS_NUMBERS
    x = X(step, :);
    y = Y(step, :);
    P_before = x' * x;
    [A, P_now, K, sigma, e_hat] = ExpWinRLS(x, y, P_before, A, lambda);
    error(i, :) = e_hat;
    i = i + 1;
end

figure();
plot(1:1:L-2*AR_PARAMETERS_NUMBERS+1, error);
title('error')
xlabel('dots');
ylabel('amplitude');
grid on;

predict_numbers = 500;
y_predict = zeros([predict_numbers, 1]);
for j=1:AR_PARAMETERS_NUMBERS
    y_predict(j, :) = Y(j,:);
end
for i= AR_PARAMETERS_NUMBERS+1:predict_numbers
    x_predict = zeros([1, AR_PARAMETERS_NUMBERS]);
    for j= 1:AR_PARAMETERS_NUMBERS
        x_predict(1, j) = y_predict(i-j);
    end
    y_predict(i, :) = x_predict * A;
end
figure();
plot((1:1:predict_numbers)*Ts, y_predict);
title('simulate model estimate signal')
xlabel('time (t)');
ylabel('amplitude');
grid on;

A

% ar1 = 2cos(2*w*delta_t)
freq_id = acos(A(1)/2) / (delta_t*2*pi)


% state space Ad
Ad = [1 0;A']
eigvalue = eig(Ad)



