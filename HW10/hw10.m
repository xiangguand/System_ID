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
t_begin = 10;
lambda = 0.95;  % forgeting factor
t_size = t_begin;
% input, rand generate
X = rand([L, p]);
X(:, 1) = X(:, 1) * 2 + 1;
X(:, 2) = X(:, 2) / 5 + 3.2;
X(:, 3) = X(:, 3) * 3.14 - 0.5;
% X(:, 4) = X(:, 3)+0.001;
X(:, 4) = X(:, 4) * 90 + 3.25;

% As we know A matrix, this is the answer
A = [1 2 -5;9 1.21 5;7 -2 3.14;5 2 -2.02;8 7 6.58]    % pxq

% random noise
E_n = rand([L, q])*10;

% output
Y = X * A + E_n;

%% step 3, use LSE
A_y_useLSE = pseudo_inverse(X'*X)*X'*Y
E_y_useLSE = Y - X * A_y_useLSE;

%% step 2
A_x = rand([p, p]);
A_y = rand([p, q]);
W = A_y;
erRLS_array = zeros([L-t_size+1, 1]);
erQR_array = zeros([L-t_size+1, 1]);
erLSE_array = zeros([L-t_size+1, 1]);
c = 1;
xx = 1:1:t_size;
% using RLS (Recusive Least Square) to tune our parameters
for t = t_begin+1:L
    % channel initial
    e_x = X(t-t_begin:t-1, :);   % n_samples x Dimension_input
    e_y = Y(t-t_begin:t-1, :);   % n_samples x Dimension_output
    d = e_y;
    [Q, R] = GramSchmidtQR(e_x);
    phi = e_x'*(lambda*eye(t_size))*e_x;
    yRLS = e_x*A_y;
    erRLS = lambda*(d - yRLS);
    A_y = A_y + pseudo_inverse(phi)*e_x'*erRLS;
    erRLS_array(c, 1) = sum(sum(abs(erRLS)))/t_size;
    
    U = Q'*e_x;
    w = pseudo_inverse(U)*Q'*e_y;
    erQR = lambda*(d - e_x*w);
    erQR_array(c, 1) = sum(sum(abs(erQR)))/t_size;
    yQR = e_x*w;
    
    yLSE = e_x*A_y_useLSE;
    erLSE = d - yLSE;
    erLSE_array(c, 1) = sum(sum(abs(erLSE)))/t_size;
    c = c + 1;

end

xdot = linspace(1, L-t_size+1, L-t_size+1);
figure();
hold on; grid on;
plot(xdot, erRLS_array, 'r');
plot(xdot, erLSE_array, 'b');
plot(xdot, erQR_array, 'g');
ylim([0, 15]);
legend('RLS', 'LSE', 'QRD');
xlabel('x dots');
ylabel('error');
title('error record');


total_erRLS = sum(erRLS_array)
total_erLSE = sum(erLSE_array)
total_erQR = sum(erQR_array)

%% step 4


%% step 5
