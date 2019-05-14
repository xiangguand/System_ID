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

L = 30000;    % Create 2000 terms of markov parameters

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
U = zeros([m*L, L]);
save_U = zeros([2,L]);
random_level = 2;
step = 1;
i = 1;
while step < m*L
   if step == 1
       Y(:, step:step+m-1) = Dd;
   else
       Y(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
   end
   U(step:step+m-1, i) = rand([2, 1])*random_level;  % random input
   save_U(:, i) = U(step:step+m-1, i);
   i = i + 1;
   step = step + m;
end

% response of 2 output sequences.
% change input random level, responses have the same effect.
y = Y*U;
t_dot = linspace(0, L*Ts, L);

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
chop_U = reshape(save_U, [m, L/chop_sections, chop_sections]);
chop_y = reshape(y, [2, L/chop_sections, chop_sections]);

% Do circular correlation algorithm, R_yu and R_uu in each section
N = double(int32(L/chop_sections/2));
chop_Ryu = zeros([m, N, chop_sections]);
chop_Ruu = zeros([q, N, chop_sections]);
for chop_step = 1:chop_sections
    for i=1:N
        temp_yu = zeros([q, 1]);
        temp_uu = zeros([m, 1]);
        for j=1:N
            temp_yu = temp_yu + chop_y(:, j+i, chop_step) .* chop_U(:, j, chop_step);
            temp_uu = temp_uu + chop_U(:, j+i, chop_step) .* chop_U(:, j, chop_step);
        end
        chop_Ryu(:, i, chop_step) = temp_yu / N;
        chop_Ruu(:, i, chop_step) = temp_uu / N;
    end
end
chop_Syu = zeros([m, N, chop_sections]);
chop_Suu = zeros([q, N, chop_sections]);
for i=1:chop_sections
    % input one
    chop_Syu(1,:,i) = abs(fft(chop_Ryu(1,:,i)));
    chop_Suu(1,:,i) = abs(fft(chop_Ruu(1,:,i)));
    % input two
    chop_Syu(2,:,i) = abs(fft(chop_Ryu(2,:,i)));
    chop_Suu(2,:,i) = abs(fft(chop_Ruu(2,:,i)));
end
mean_Syu = mean(chop_Syu(:, 1:N/2, :), 3);
mean_Suu = mean(chop_Suu(:, 1:N/2, :), 3);
freq_x = linspace(0, Fs/2, N/2);
figure();
subplot(2,1,1);
plot(freq_x, mean_Syu(1, :));
title('Mean of density spectrum S_{yu}1');
grid on;
subplot(2,1,2);
plot(freq_x, mean_Syu(2, :));
title('Mean of density spectrum S_{yu}2');
grid on;

figure();
subplot(2,1,1);
plot(freq_x, mean_Suu(1, :));
title('Mean of density spectrum S_{uu}1');
grid on;
subplot(2,1,2);
plot(freq_x, mean_Suu(2, :));
title('Mean of density spectrum S_{uu}2');
grid on;


% step 4, calculate the transfer function and the FIR weighting sequence
Gz = abs(mean_Syu) .* (abs(mean_Suu).^-1);
figure();
subplot(2,1,1);
plot(freq_x, Gz(1, :));
title('Mean of density spectrum G_z1');
grid on;
subplot(2,1,2);
plot(freq_x, Gz(2, :));
title('Mean of density spectrum G_z2');
grid on;

% step 5, use ifft of function of Matlab to transform the elements of G[n]
% matrix back the Markov parameters.
new_Y = [ifft(Gz(1, :)); ifft(Gz(2, :))];
new_t = linspace(0, Ts*L/chop_sections, N/2);
figure();
subplot(2,1,1);
plot(new_t, new_Y(1, :));
title('New Markov parameters Y^*_1');
grid on;
subplot(2,1,2);
plot(new_t, new_Y(2, :));
title('New Markov parameters Y^*_2');
grid on;

% % Do fft
% fft_chop_U = zeros(m, L/chop_sections, chop_sections);
% fft_chop_y = zeros(2, L/chop_sections, chop_sections);
% temp_l = 1;
% x = linspace(0, Ts*100, 100);
% for i=1:chop_sections
%     % input one
%     fft_chop_U(1,:,i) = abs(fft(chop_U(1,:,i)));
%     fft_chop_y(1,:,i) = abs(fft(chop_y(1,:,i)));
%     % input two
%     fft_chop_U(2,:,i) = abs(fft(chop_U(2,:,i)));
%     fft_chop_y(2,:,i) = abs(fft(chop_y(2,:,i)));
% end
% % mean the input and output, after fft get half of data that it is valid.
% mean_fft_U = mean(fft_chop_U(:, 1:L/chop_sections/2, :), 3);
% mean_fft_y = mean(fft_chop_y(:, 1:L/chop_sections/2, :), 3);

