clc, clear, close all
%{
    System ID Homework 10
    input dimension p = n = 5
    output dimension q = 3
%}

%% step 1
p = 5;    % input
n = p;
q = 3;    % output
L = 200;
H_crop_dim = 2;
t_begin = p;
t_size = t_begin;
% input, rand generate
X = rand([L, p]);
X(:, 1) = X(:, 1) * 2 + 1;
X(:, 2) = X(:, 2) / 5 + 3.2;
X(:, 3) = X(:, 3) * 3.14 - 0.5;
% X(:, 4) = X(:, 4) * 9 + 3.25;
% 中間一段接近singular
X(0.25*L:0.5*L, 4) = X(0.25*L:0.5*L, 3) + 0.0001;

% As we know A matrix, this is the answer
A = [1 2 -5;9 1.21 5;7 -2 3.14;5 2 -2.02;8 7 6.58]    % pxq

% random noise
E_n = rand([L, q]);
X = X';
% output
d = X'*A + E_n;

% using least spuare, it can not inverse the Rxx matrix, because it's
% singular
A_use_LS = (X*X')^-1*X*d

lambda = 0.9;   % forgeting factor
epsilon = 0.1;
inv_R = epsilon*eye(p);
w = randn([p, q]);
er = zeros([L, 2]);
y_pred = zeros([L, q]);
LS_pred = zeros([L, q]);
for i=1:L
    e = d(i, :) - X(:, i)' * w;
    e2 = d(i, :) - X(:, i)' * A_use_LS;
    y_pred(i, :) = X(:, i)' * w;
    LS_pred(i, :) = X(:, i)' * A_use_LS;
    k = inv_R*X(:, i);
    kappa = k / (lambda + X(:, i)'*k);
    inv_R = 1/lambda*(inv_R - ((kappa*kappa')/(lambda + X(:, i)'*k)));
    w = w + kappa*e;
    % collect error
    er(i, 1) = mean(e.^2);
    er(i, 2) = mean(e2.^2);
end
w
er_sum = sum(er)

figure();
plot(er(:, 1));
hold on; grid on;
plot(er(:, 2));
legend('RLS', 'LS');
title('Loss');
xlabel('dots');
ylabel('MSE');

figure();
plot(y_pred(:, 1), 'r^--');
hold on; grid on;
plot(LS_pred(:, 1), 'b>--');
plot(d(:, 1), 'ko--');
legend('RLSpredict', 'LSpredict', 'true');
title('Prediction & True value (output 1)');
xlabel('dots');
ylabel('value');

figure();
plot(y_pred(:, 2), 'r^--');
hold on; grid on;
plot(LS_pred(:, 2), 'b>--');
plot(d(:, 2), 'ko--');
legend('RLSpredict', 'LSpredict', 'true');
title('Prediction & True value (output 2)');
xlabel('dots');
ylabel('value');

figure();
plot(y_pred(:, 3), 'r^--');
hold on; grid on;
plot(LS_pred(:, 3), 'b>--');
plot(d(:, 3), 'ko--');
legend('RLSpredict', 'LSpredict', 'true');
title('Prediction & True value (output 3)');
xlabel('dots');
ylabel('value');


