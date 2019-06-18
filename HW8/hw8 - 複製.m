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
chop_sections = 2;

L = 2000;    % Create 2000 terms of markov parameters

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
U_noise = zeros([m*L, L]);
save_noise = zeros([m,L]);
U_unit = zeros([m*L, L]);
random_level = 1;
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = Dd;
   else
       Y(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
   end
   U_noise(step:step+m-1, i) = rand([2, 1])*random_level;  % random input
   save_noise(:, i) = U_noise(step:step+m-1, i);
   U_unit(step:step+m-1, i) = u0;
   i = i + 1;
   step = step + m;
end
noice = save_noise;
% response of 2 output sequences.
% change input random level, responses have the same effect.
y = Y*(U_noise+U_unit);
y_markov = Y*U_unit;
t_dot = linspace(0, L*Ts, L);

figure();
subplot(2,1,1);
plot(t_dot, save_noise(1,:));
title("output1 noise");
xlabel("time (s)");
ylabel("Value");
grid on;

subplot(2,1,2);
plot(t_dot, save_noise(2,:));
title("output2 noise");
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

Ryu = zeros([m, L]);
Ruu = zeros([m, L]);
Ryy = zeros([q, L]);
Syu = zeros([m, L]);
Suu = zeros([m, L]);
Syy = zeros([q, L]);
for j=1:q
     temp = xcorr(y(j, :), save_noise(j, :), 'biased');
     Ryu(j,:) = temp(:, 1:L);
     Syu(j, :) = fft(Ryu(j, :));
     temp = xcorr(y(j, :), y(j, :), 'biased');
     Ryy(j,:) = temp(:, 1:L);
     Syy(j, :) = fft(Ryy(j, :));
end
for j=1:m
    temp = xcorr(save_noise(j, :), save_noise(j, :), 'biased');
    Ruu(j,:) = temp(:, 1:L);
    Suu(j, :) = fft(Ruu(j, :));
end

% plot the R_yu and R_uu
x_t = linspace(0, L*Ts, L);
figure();
subplot(2,1,1);
plot(x_t, Ryu(1, :));
title('R_{yu}1');
grid on;
subplot(2,1,2);
plot(x_t, Ryu(2, :));
title('R_{yu}2');
grid on;

figure();
subplot(2,1,1);
plot(x_t, Ruu(1, :));
title('R_{uu}1');
grid on;
subplot(2,1,2);
plot(x_t, Ruu(2, :));
title('R_{uu}2');
grid on;

figure();
subplot(2,1,1);
plot(x_t, Ryy(1, :));
title('R_{yy}1');
grid on;
subplot(2,1,2);
plot(x_t, Ryy(2, :));
title('R_{yy}2');
grid on;

% step 3, chop the input and output sequences for several sections
% Do circular correlation algorithm, R_yu and R_uu in each section
N = 1 + (L/chop_sections-1)*2;
chop_Ryu = zeros([m, L/chop_sections, chop_sections]);
chop_Ruu = zeros([q, L/chop_sections, chop_sections]);
chop_Ryy = zeros([q, L/chop_sections, chop_sections]);
chop_Syu = zeros([m, L/chop_sections, chop_sections]);
chop_Suu = zeros([q, L/chop_sections, chop_sections]);
chop_Syy = zeros([q, L/chop_sections, chop_sections]);
for chop_step = 1:chop_sections
    for j=1:q
        chop_Ryu(j, :, chop_step) = Ryu(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));
        chop_Ruu(j, :, chop_step) = Ruu(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));
        chop_Ryy(j, :, chop_step) = Ryy(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));

        chop_Syu(j,:,chop_step) = Syu(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));
        chop_Suu(j,:,chop_step) = Suu(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));
        chop_Syy(j,:,chop_step) = Syy(j, (chop_step-1)*(L/chop_sections)+1:(chop_step)*(L/chop_sections));
    end
end

mean_Syu = mean(chop_Syu(:, 1:L/chop_sections, :), 3);
mean_Suu = mean(chop_Suu(:, 1:L/chop_sections, :), 3);
mean_Syy = mean(chop_Syy(:, 1:L/chop_sections, :), 3);
freq_x = linspace(0, Fs/2, L/chop_sections);
figure();
subplot(2,1,1);
plot(freq_x, abs(mean_Suu(1, :)));
title('Mean of density spectrum S_{uu}1');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;
subplot(2,1,2);
plot(freq_x, abs(mean_Suu(2, :)));
title('Mean of density spectrum S_{uu}2');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

figure();
subplot(2,1,1);
plot(freq_x, abs(mean_Syu(1, :)));
title('Mean of density spectrum S_{yu}1');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;
subplot(2,1,2);
plot(freq_x, abs(mean_Syu(2, :)));
title('Mean of density spectrum S_{yu}2');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;


% step 4, calculate the transfer function and the FIR weighting sequence
Gz = mean_Syu .* (mean_Suu.^-1);
figure();
subplot(2,1,1);
plot(freq_x, abs(Gz(1, :)));
title('Mean of density spectrum G_z1');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;
subplot(2,1,2);
plot(freq_x, abs(Gz(2, :)));
title('Mean of density spectrum G_z2');
xlabel("frequency (Hz)");
ylabel("Amplitude");
grid on;

% step 5, use ifft of function of Matlab to transform the elements of G[n]
% matrix back the Markov parameters.
new_Y = [ifft(Gz(1, :)); ifft(Gz(2, :))];
new_t = linspace(0, Ts*L/chop_sections, L/chop_sections);
figure();
subplot(2,1,1);
plot(new_t, abs(new_Y(1, :)));
title('New Markov parameters Y^*_1');
xlabel('seconds');
grid on;
subplot(2,1,2);
plot(new_t, abs(new_Y(2, :)));
title('New Markov parameters Y^*_2');
xlabel('seconds');
grid on;

% step 6, plot the Bode plot of spectrum sequences and Mean Spectrum
% sequence and discuss to discuss the effect of SNR
SNR1 = abs(mean_Syu(1,:))/(abs(mean_Syu(1,:)./mean_Suu(1,:)))
SNR2 = abs(mean_Syu(2,:))/(abs(mean_Syu(2,:)./mean_Suu(2,:)))


