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
L = 20;
H_crop_dim = 2;
t_begin = p;
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
N = p;
delta = 0.1;
U = fliplr(eye(p)*delta);
w = zeros([p, q]);
dq2 = zeros([1, q]);
mygamma = 1;
for k=2:L
    xaux = X(k, :);
    Uaux = lambda^0.5*U;
    daux = Y(k, :)';     % qx1
    dq2aux = lambda^0.5*dq2;
    gamma = 1;
    for n=1:N
        cosTheta = abs(U(N-n+1, n)) / sqrt(xaux(n)^2 + Uaux(N-n+1, n)^2);
        sinTheta = xaux(n) / Uaux(N-n+1, n)*cosTheta;
        xaux(n) = 0;
        Uaux(N-n+1, n) = sinTheta*xaux(n) + cosTheta*Uaux(N-n+1, n);
        for m=n+1:N 
            oldxaux = xaux(m);
            xaux(m) = cosTheta*oldxaux - sinTheta*Uaux(N-n+1, m);
            Uaux(N-n+1, m) = sinTheta*oldxaux + cosTheta*Uaux(N-n+1, m)
        end
        % Obtain \gamma
        gamma = gamma*cosTheta;
        % Obtain eq1(k) and updating dq2(k)
        olddaux = daux;
        daux = cosTheta*olddaux - sinTheta*dq2aux;
        dq2aux = sinTheta*olddaux - sinTheta*dq2aux
        % Back-substitution, if necessary
        summa = 0;
        for m=1:n-2
            summa = summa + Uaux(n, N-m+1);
        end
        dq2aux
        w = (dq2aux(n) - summa) / Uaux(n , N-n+1)
    end
    U = Uaux;
    mygamma = gamma;
    dq2 = dq2aux;
    eq1 = daux;
    epsilon = eq1'*mygamma;
    e = eq1'*mygamma;
end



%% step 4


%% step 5
