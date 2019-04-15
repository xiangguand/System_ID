clc, clear, close all
%{
  @System ID Homework5
  @Please write a program which can transform the representation of a
  system from state space model to the linear defference model.
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 100
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    A_bar matrix: nxn
    B matrix: nxm
    B_bar matrix: nx(m+q)
    G matrix: nxq
    C matrix: qxn
    D matrix: qxm
   
    y = h * V    % h is markov parameters
    y: qxL
    h: qx[(q+m)(L-1)+m] = Y (markov parameters) or call it weighting
       sequence.
    V: [(q+m)(L-1)+m]xL
    v: (m+q)x1

    The linear difference model is also called ARX model (Auto-Regressing
    eXogenous model)
    % The final equation is showed below:
    y[k] = sum(Yi2*y[k-i], i=1->4) + sum(Yi1*u[k-i], i=1->4) + D*u[k]
%}

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
% Define dimension
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;
delta_t = 0.01;
SYSTEM_DIM = 4;
L = 100;    % Create 100 terms of markov parameters
u0 = [1;1];
y0 = [0;0];
v0 = [u0;y0];

% Design observer feedback gain
Gc = [0 0;0 0;0 0;0 0];
Gd = [0 0;0 0;0 0;0 0];

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t)

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

Gd(:,1) = -Ad(:,1);
Gd(:,2) = -Ad(:,3);
Ad_bar = Ad + Gd*Cd
Bd_bar = [Bd + Gd*Dd -Gd]
eig_Ad_bar = eig(Ad_bar)

% using ARX model
Yi_1 = zeros([q, m*SYSTEM_DIM]);
Yi_2 = zeros([q, m*SYSTEM_DIM]);
u = zeros([m, L]);
u(:, 1) = u0;
y = zeros([m, L]);
y(:, 1) = y0;

% i=1 ==> markov parameters is equal to D. Dimension is different.
step = 1;
i = 1;
while step < m*SYSTEM_DIM
    if step==1
        temp_data = [Dd [0 0;0 0]];
    else
        temp_data = Cd * Ad_bar^(i-1) * Bd_bar;    % 2x4
    end
    Yi_1(:, step:step+m-1) = temp_data(:,1:2); % 2x2
    Yi_2(:, step:step+m-1) = temp_data(:,3:4); % 2x2
    step = step + 2;
    i = i + 1;
end
% Do ARX.
for k=5:L
    temp_data2 = zeros([m, 1]);
    temp_data1 = zeros([m, 1]);
    s = 1;
    for i=1:4
        temp_data2  = temp_data2 + Yi_2(:,s:s+1)*y(:, k-i);
        temp_data1 =  temp_data1 + Yi_1(:,s:s+1)*u(:, k-i);
        s = s + 2;
    end
    y(:, k) = temp_data2 + temp_data1 + Dd*u(:,k);
end

x_dot = linspace(1, L, L);
figure();
plot(x_dot, y(1, :));
xlabel('k samples');
ylabel('output');
title('y1 Markov parameters by using ARX model');
grid on;

figure();
plot(x_dot, y(2, :));
xlabel('k samples');
ylabel('output');
title('y2 Markov parameters by using ARX model');
grid on;

