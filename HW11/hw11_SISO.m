clc, clear, close all

u0 = 1;
m = 1;
q = 1;
n = 2;
L = 1000;
alpha = 5;

A = [0 1;1 1];
B = [1;1];
C = [1 0];
D = [0];
[num, den] = ss2tf(A, B, C, D);
T_original = tf(num, den)

% Calculate system Markov parameters
% using symbol hn(k)
Y = zeros([q, m*L]);
U = zeros([m*L, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = D;
   else
       Y(:, step:step+m-1) = C*(A^(i-2))*B;
   end
   U(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end

yn = Y*U;

% Constuct the Hankel matrix H0
H0 = zeros([alpha, alpha]);
% オ场
for i=1:alpha
    Hn = Y(:, i+1);
    for j=1:i
        index_r = j;
        index_c = (i-j+1);
        H0(index_r, index_c) = Hn;
    end
end
% 场
for i=m+1:alpha
    Hn = Y(:, alpha+i);
    for j=alpha:-1:i
        index_r = i+alpha-j;
        index_c = j;
        H0(index_r, index_c) = Hn;
    end
end

% Constuct the Hankel matrix H0
H1 = zeros([alpha, alpha]);
% オ场
for i=1:alpha
    Hn = Y(:, i+2);
    for j=1:i
        index_r = j;
        index_c = (i-j+1);
        H1(index_r, index_c) = Hn;
    end
end
% 场
for i=m+1:alpha
    Hn = Y(:, alpha+i+1);
    for j=alpha:-1:i
        index_r = i+alpha-j;
        index_c = j;
        [index_r index_c]
        H1(index_r, index_c) = Hn;
    end
end
% EVD the Markov parameters
[R, Sigma, S] = svd(H0);
O = R*sqrtm(Sigma)      % observable matrix
C = sqrtm(Sigma)*S'     % controllable matrix
A_bar = Sigma^-0.5*R'*H1*S*Sigma^-0.5

A_hat = A_bar(1:2, 1:2)
B_hat = C(1:2, 1)
C_hat = O(1, 1:2)
D = [0];

[num, den] = ss2tf(A_hat, B_hat, C_hat, D);
T_svd = tf(num, den)


% T_orginal should the same as T_svd

T_original
T_svd

pole_original = pole(T_original);
pole_svd = pole(T_svd);
figure();
plot(pole_original, 'b*');
hold on;grid on;
plot(pole_svd, 'ro');
legend('original', 'EVD');




