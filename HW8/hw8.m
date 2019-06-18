    
clc, clear, close all
%{
  @System ID Homework8
  @Use pseudo random generator to generate excitation input sequences for
  the two forcing input
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 100
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    B matrix: nxm
    C matrix: qxn
    D matrix: qxm
   
    y = Y * U    % h is markov parameters
    Y: qxm
    u: mxL
%}
% Define dimension
Fs = 10;
Ts = 1/Fs;
delta_t = Ts;
chop_sections = 5;

L = 1000;    % Create 1000 terms of markov parameters

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
u0 = para_struct.input_vect;
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;

ssd = ss(Ac, Bc, Cc, Dc);
sysd = c2d(ssd, delta_t);

Ad = expm(Ac*delta_t);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;

e = eig(Ad)

% Calculate system Markov parameters (symbol h)
% using symbol hn(k)
Y = zeros([q, m*L]);
% Random noise generates, to exite system
noise = zeros([m*L, L]);
U = zeros([m*L, L]);
save_noise = zeros([2,L]);
save_U = zeros([2,L]);
random_level = 20;
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = Dd;
   else
       Y(:, step:step+m-1) = Cd*(Ad^(i-2))*Bd;
   end
   noise(step:step+m-1, i) = rand([2, 1])*random_level;  % random input
   U(step:step+m-1, i) = u0;
   save_noise(:, i) = noise(step:step+m-1, i);
   save_U(:, i) = u0;
   i = i + 1;
   step = step + m;
end

% response of 2 output sequences.
% change input random level, responses have the same effect.
y = Y*noise;
t_dot = linspace(0, L*Ts, L);

figure();
subplot(2,1,1);
plot(t_dot, save_noise(1,:));
title("Random input data 1");
xlabel("time (s)");
ylabel("Value");
grid on;

subplot(2,1,2);
plot(t_dot, save_noise(2,:));
title("Random input data 2");
xlabel("time (s)");
ylabel("Value");
grid on;

figure();
subplot(2,1,1);
plot(t_dot, y(1,:));
title("output1 random exitation");
xlabel("time (s)");
ylabel("Value");
grid on;

subplot(2,1,2);
plot(t_dot, y(2,:));
title("output2 random exitation");
xlabel("time (s)");
ylabel("Value");
grid on;

y1_fft = abs(fft(y(1, :)));
y2_fft = abs(fft(y(2, :)));
y1_fft = y1_fft(1:L/2);
y2_fft = y2_fft(1:L/2);

freq_dot = linspace(0, Fs/2, L/2);
figure();
subplot(2,1,1);
plot(freq_dot, y1_fft);
title("output1 random exitation");
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

subplot(2,1,2);
plot(freq_dot, y2_fft);
title("output2 random exitation");
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

% step 3, chop the input and output sequences for several sections
chop_Y = reshape(Y, [2, 2*L/chop_sections, chop_sections]);
chop_U = reshape(save_noise, [m, L/chop_sections, chop_sections]);
chop_y = reshape(y, [q, L/chop_sections, chop_sections]);

% Do circular correlation algorithm, R_yu and R_uu in each section
N = 1 + (L/chop_sections-1)*2;
chop_Ryu = zeros([m, N, chop_sections]);
chop_Ruu = zeros([q, N, chop_sections]);
chop_Sy = zeros([m, L/chop_sections, chop_sections]);
chop_Su = zeros([q, L/chop_sections, chop_sections]);
for chop_step = 1:chop_sections
    for j=1:q
        chop_Ryu(j, :, chop_step) = xcorr(chop_y(j, :, chop_step), chop_U(j, :, chop_step), 'biased');
        chop_Ruu(j, :, chop_step) = xcorr(chop_U(j, :, chop_step), chop_U(j, :, chop_step), 'biased');
        
%         chop_Syu(j,:,chop_step) = fft(chop_Ryu(j,:,chop_step));
%         chop_Suu(j,:,chop_step) = fft(chop_Ruu(j,:,chop_step));
        chop_Sy(j,:,chop_step) = fft(chop_y(j, :, chop_step));
        chop_Su(j,:,chop_step) = fft(chop_U(j, :, chop_step));
    end
end

% plot the R_yu and R_uu
flat_Ryu = reshape(chop_Ryu, [q, N*chop_sections]);
flat_Ruu = reshape(chop_Ruu, [q, N*chop_sections]);
flat_Sy = reshape(chop_Sy, [q, L]);
flat_Su = reshape(chop_Su, [q, L]);
x_t = linspace(0, N*chop_sections*Ts, N*chop_sections);
figure();
subplot(2,1,1);
plot(x_t, abs(flat_Ryu(1, :)));
title('flat R_{yu}1');
grid on;
subplot(2,1,2);
plot(x_t, abs(flat_Ryu(2, :)));
title('flat R_{yu}2');
grid on;

figure();
subplot(2,1,1);
plot(x_t, abs(flat_Ruu(1, :)));
title('flat R_{uu}1');
grid on;
subplot(2,1,2);
plot(x_t, abs(flat_Ruu(2, :)));
title('flat R_{uu}2');
grid on;

figure();
subplot(2,1,1);
plot(abs(flat_Sy(1, :)));
title('flat S_{y}1');
grid on;
subplot(2,1,2);
plot(abs(flat_Sy(2, :)));
title('flat S_{y}2');
grid on;

figure();
subplot(2,1,1);
plot(abs(flat_Su(1, :)));
title('flat S_{u}1');
grid on;
subplot(2,1,2);
plot(abs(flat_Su(2, :)));
title('flat S_{u}2');
grid on;

mean_Ryu = mean(chop_Ryu(:, 1:N, :), 3);
mean_Ruu = mean(chop_Ruu(:, 1:N, :), 3);
x_t = linspace(0, N*Ts, N);
figure();
subplot(2,1,1);
plot(x_t, mean_Ruu(1, :));
title('Mean of R_{uu}1');
xlabel("time (s)");
ylabel("Value");
grid on;
subplot(2,1,2);
plot(x_t, mean_Ruu(2, :));
title('Mean of R_{uu}2');
xlabel("time (s)");
ylabel("Value");
grid on;

figure();
subplot(2,1,1);
plot(x_t, mean_Ryu(1, :));
title('Mean of R_{yu}1');
xlabel("time (s)");
ylabel("Value");
grid on;
subplot(2,1,2);
plot(x_t, mean_Ryu(2, :));
title('Mean of R_{yu}2');
xlabel("time (s)");
ylabel("Value");
grid on;


mean_Sy = mean(chop_Sy(:, :, :), 3);
mean_Su = mean(chop_Su(:, :, :), 3);
freq_x = linspace(0, Fs/2, L/chop_sections);
figure();
subplot(2,1,1);
plot(freq_x, abs(mean_Su(1, :)));
title('Mean of density spectrum S_{uu}1');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;
subplot(2,1,2);
plot(freq_x, abs(mean_Su(2, :)));
title('Mean of density spectrum S_{uu}2');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

figure();
subplot(2,1,1);
plot(freq_x, abs(mean_Sy(1, :)));
title('Mean of density spectrum S_{yu}1');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;
subplot(2,1,2);
plot(freq_x, abs(mean_Sy(2, :)));
title('Mean of density spectrum S_{yu}2');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;


% step 4, calculate the transfer function and the FIR weighting sequence
Gz = mean_Sy .* (mean_Su.^-1);
figure();
subplot(2,1,1);
plot(freq_x, abs(Gz(1, :)));
title('Mean of density spectrum G_z1');
grid on;
subplot(2,1,2);
plot(freq_x, abs(Gz(2, :)));
title('Mean of density spectrum G_z2');
grid on;

% step 5, use ifft of function of Matlab to transform the elements of G[n]
% matrix back the Markov parameters.
new_Y = [ifft(Gz(1, :)); ifft(Gz(2, :))];
new_t = linspace(0, Ts*L/chop_sections, L/chop_sections);
figure();
subplot(2,1,1);
plot(new_t, abs(new_Y(1, :)));
title('New Markov parameters Y^*_1');
grid on;
subplot(2,1,2);
plot(new_t, abs(new_Y(2, :)));
title('New Markov parameters Y^*_2');
grid on;