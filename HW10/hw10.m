clc, clear, close all
%{
    System ID Homework 10
    input dimension p = 5
    output dimension q = 2
    input matrix H
    output matrix Y
    ID parameters A_y
    X: txp
    Y: txq
    H: txchop_num
d
    Hi: txi 
    Ei_y: txq
    Ei_x: tx(p-i)
    Ai_y: i x q
    Ai_x: i x p
    A_y: p x q
    A_x: p x p
    data length L
    Y^{hat} = H*A_y
    E_x = X - H*A_x
    E_y = Y - H*A_y
%}

%% step 1
p = 5;
q = 2;
L = 2000;
H_crop_dim = 2;
t_begin = 6;
lambda = 1;
t_size = t_begin;
% input, rand generate
X = rand([L, p]);
X(:, 1) = X(:, 1) * 2;
X(:, 2) = X(:, 2) / 5;
X(:, 3) = X(:, 3) * 3.14;
X(:, 4) = X(:, 4) * 90;

H = X(:, 1:H_crop_dim);

% As we know A matrix, this is the answer
A = [1 2;9 1.21;7 -2;5 2;8 7]

% random noise
E_n = rand([L, q])*0.01;

% output
Y = X * A + E_n;

%% step 2
% initial
A_x = rand([p, p, p]);
A_y = rand([p, q, p]);

R_xx = zeros([p, p]);
R2 = zeros([p, p]);
R_xy = zeros([q, p]);
R_uu = pining(p);
R3 = R_uu;
C_xx = zeros([p, p]);
C_xy = zeros([p, q]);
% using RLS (Recusive Least Square) to tune our parameters
for t = t_begin:L
    % channel initial
    e_x = X(t-t_size+1:t, :)';   % X dim x n_samples 
    e_x2 = e_x
    e_y = Y(t-t_size+1:t, :)';   % Y dim x n_samples
    for i=2:p   % channel update
        R_uu(i) = R_uu(i-1) - e_x(i, i-1)*pseudo_inverse(R_xx(i, i-1))*e_x(i, i-1);
        R_xx(:, i-1) = lambda*R_xx(:, i-1) + e_x(i, i-1)*pseudo_inverse(R_uu(i-1))*e_x(:, i-1);
        e_x(:, i) = e_x(:, i-1) + pseudo_inverse(R_xx(i, i-1))*R_xx(:, i-1);
        R_xy(:, i-1) = lambda*R_xy(:, i-1) + e_x(i, i-1)*pseudo_inverse(R_uu(i-1))*e_y(:, i-1);
        e_y(:, i) = e_y(:, i) - e_x(i, i-1)*pseudo_inverse(R_xx(i, i-1))*R_xy(:, i-1);
        R3(i) = R3(i-1) - e_x2(i, i-1)*pseudo_inverse(R2(i, i-1))*e_x2(i, i-1);
        for j=i+1:p
            R2(j, i-1) = lambda*R2(j, i-1) + e_x2(i, i-1)*pseudo_inverse(R3(i-1))*e_x2(j, i-1);
            e_x2(j, i) = e_x2(j, i-1) - e_x2(j, i-1)*pseudo_inverse(R2(i, i-1))*R2(j, i-1);
        end
    end % end i
    for i=2:p
        for j=i+1:p
            C_xx(j, i-1) = pseudo_inverse(R_xx(i, i-1))*R_xx(j, i-1); 
            temp_Ax = A_x(:, :, i);
            temp_Ax(end, :) = 0;
            temp2_Ax = A_x(:, i, i-1);
            temp2_Ax(end, :) = -1;
            A_x(:, :, i) = temp_Ax - temp2_Ax*C_xx(i-1, :);
        end
    end
    for i=2:q
        for j=1:q
            C_xy(j, i-1) = pseudo_inverse(R_xx(i, i-1))*R_xy(j, i-1);
            temp_Ay = A_y(:, :, i);
            temp_Ay(end, :) = 0;
            temp2_Ay = A_y(:, i, i-1);
            temp2_Ay(end, :) = -1;
            A_y(:, :, i) = temp_Ay - temp2_Ay*C_xy(i-1, :);
        end
    end
    y_hat = Y(t, :)' - e_y(:, i);
    err_y = e_y' - e_x'*A_y(:,:,end);
    err_x = e_x' - e_x'*A_x(:,:,end);

end


%% step 3, use LSE
A_y_useLSE = (X'*X)^-1*X'*Y
E_y_useLSE = Y - X * A_y_useLSE;


%% step 4


%% step 5

