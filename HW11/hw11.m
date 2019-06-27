clc, clear, close all

% Define dimension
Fs = 10;
Ts = 1/Fs;
delta_t = Ts;
chop_sections = 5;
alpha = 4;

L = 1000;    % Create 1000 terms of markov parameters

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
u0 = para_struct.input_vect;
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t);
T_original = tf(sysd);

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

e = eig(Ad)

% Calculate system Markov parameters
% using symbol hn(k)
Y = zeros([q, m*L]);
U = zeros([m*L, L]);
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = Dd;
   else
       Y(:, step:step+m-1) = Cd*(Ad^(i-2))*Bd;
   end
   U(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end

yn = Y*U;

% Constuct the Hankel matrix H0
H0 = zeros([alpha*m, alpha*q]);

% オ场
for i=1:m:alpha*m
    Hn = Y(:, i*(m-1)+m:i*(m-1)+m+1);
    [i*(m-1)+m, i*(m-1)+m+1];
    for j=1:q:i
        index_r = j;
        index_c = (i-j+1);
        H0(index_r:index_r+1, index_c:index_c+1) = Hn;
    end
end
% 场
for i=m+1:m:alpha*m
    Hn = Y(:, alpha+i*(m-1)+2*m:alpha+i*(m-1)+2*m+1);
    [alpha+i*(m-1)+2*m, alpha+i*(m-1)+2*m+1];
    for j=alpha*m-1:-q:i
        index_r = i+(alpha*m-j)-1;
        index_c = j;
        H0(index_r:index_r+1, index_c:index_c+1) = Hn;
    end
end

% Constuct the Hankel matrix H1
H1 = zeros([alpha*m, alpha*q]);
n1 = 0;
% オ场
for i=1:m:alpha*m
    Hn = Y(:, i*(m-1)+2*m:i*(m-1)+2*m+1);
    [i*(m-1)+2*m, i*(m-1)+2*m+1];
    for j=1:q:i
        index_r = j;
        index_c = (i-j+1);
        H1(index_r:index_r+1, index_c:index_c+1) = Hn;
    end
end
% 场
for i=m+1:m:alpha*m
    Hn = Y(:, alpha+i*(m-1)+3*m:alpha+i*(m-1)+3*m+1);
    [alpha+i*(m-1)+3*m,alpha+i*(m-1)+3*m+1];
    for j=alpha*m-1:-q:i
        index_r = i+(alpha*m-j)-1;
        index_c = j;
        H1(index_r:index_r+1, index_c:index_c+1) = Hn;
    end
end
% EVD the Markov parameters
[R, Sigma, S] = svd(H0);
P = R*sqrtm(Sigma)      % observable matrix
Q = sqrtm(Sigma)*S'     % controllable matrix
A_bar = Sigma^-0.5*R'*H1*S*Sigma^-0.5;
C_hat = P(1:2, 1:4)
B_hat = Q(1:4, 1:2)
A_hat = A_bar(1:4, 1:4)
D = [0 0;0 0];
sys_EVD = ss(A_hat, B_hat, C_hat, D, 0.0000000001);
impulse(sys_EVD);
T_svd = tf(sys_EVD);

T_original
T_svd

pole_original = pole(T_original);
pole_EVD = pole(T_svd);
figure();
plot(pole_original, 'b*');
hold on;grid on;
plot(pole_EVD, 'ro');
legend('original', 'EVD');









