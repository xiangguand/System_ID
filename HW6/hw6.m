clc, clear, close all
%{
  @System ID Homework6
  @Please write a program which can transform the representation of a
  system from state space model to the linear defference model.
  @Dimensionss:
    input: m = 2
    output: q = 2
    state: n = 4
    sample: L = 1000
    input matrix: mx1
    output matrix: qx1
    A marix: nxn
    B matrix: nxm
    C matrix: qxn
    D matrix: qxm
    
    u: mxL
    u[k] = sin(w*k*delta_t)
    Yi: qxm
    y: qxL    

%}

[Ac, Bc, Cc, Dc, para_struct] = createMassDampingSpringModel('hw3_machanic_n2.json')
% Define dimension
m = para_struct.input_sz;
q = para_struct.output_sz;
n = para_struct.state_sz;
% step 1 : Decide the sampling rate according to the analog system
Fs = 100;
Ts = 1/Fs;

SYSTEM_DIM = 4;
L = 200;    % Create 200 terms of markov parameters
freq_L = 0.5;
freq_H = 40;
freq_size = 10;
time_value = (freq_H/freq_L)^(1/freq_size);

% step 2 : Degitalize the system with the above sampling rate to obtain
% its digital model
Ad = expm(Ac*Ts);
Bd = Ac^-1*(Ad-eye(4))*Bc;
Cd = Cc;
Dd = Dc;
ssd = ss(Ad, Bd, Cd, Dd);
figure();
bode(ssd);
grid on;

% step 3 : change input sinusoid wave frequency from low to high.


Yi = zeros([q, m*L]);
step = 1;
while step < m*L
   if step == 1
       Yi(:, step:step+m-1) = Dd;
   else
       Yi(:, step:step+m-1) = Cd*(Ad^(step-2))*Bd;
   end
   step = step + m;
end

% ##################################################### %
% Input 1 as sinwave, Input 2 all zero
% Becareful, sine wave frequency should be less than Fs/2
collect_Am1 = zeros([1, freq_size]);
collect_Ph1 = zeros([1, freq_size]);
collect_Am2 = zeros([1, freq_size]);
collect_Ph2 = zeros([1, freq_size]);
cnt = 1;
sin_freq = freq_L;
while true
    w = 2*pi*sin_freq;
    A = 1;
    dot1 = linspace(1, L*2, L*2);
    dot2 = zeros([1, L*2]);
    t_dot = [dot1 ; dot2];

%     figure();
%     plot(dot1, sin(w*dot1*Ts));
%     title("input");

    y = zeros([q, L]);
    for k=1:L
        step = 1;
        for tau=1:L
           u = sin(w*t_dot*Ts);
           y(:, k) = y(:, k) + Yi(:, step:step+m-1) * u(:, k + tau); 
           step = step + m;
        end
    end
    sin_freq = sin_freq * time_value;
    
    [Am1, Ph1] = getSineAmplitudeAndPhase(y(1, :));
    [Am2, Ph2] = getSineAmplitudeAndPhase(y(2, :));
    collect_Am1(1, cnt) = 20*log10(Am1);
    collect_Ph1(1, cnt) = Ph1;
    collect_Am2(1, cnt) = 20*log10(Am2);
    collect_Ph2(1, cnt) = Ph2;
    cnt = cnt + 1;
    if cnt > freq_size
        break
    end
end


freq_x = linspace(freq_L, freq_H, freq_size);

figure();
subplot(2,1,1)
plot(freq_x, collect_Am1);
title("AM1 by input1");
xlabel("Frequency (Hz)");
ylabel("db");
grid on;
set(gca, 'XScale', 'log');
subplot(2,1,2)
plot(freq_x, collect_Ph1);
title("Ph1 by input1");
xlabel("Frequency (Hz)");
ylabel("度");
grid on;
set(gca, 'XScale', 'log');

figure();
subplot(2,1,1)
plot(freq_x, collect_Am2);
title("AM2 by input1");
xlabel("Frequency (Hz)");
ylabel("db");
grid on;
set(gca, 'XScale', 'log');
subplot(2,1,2)
plot(freq_x, collect_Ph2);
title("Ph2 by input1");
xlabel("Frequency (Hz)");
ylabel("度");
grid on;
set(gca, 'XScale', 'log');

% ##################################################### %
% Input 2 as sinwave, Input 1 all zero
collect_Am1 = zeros([1, freq_size]);
collect_Ph1 = zeros([1, freq_size]);
collect_Am2 = zeros([1, freq_size]);
collect_Ph2 = zeros([1, freq_size]);
cnt = 1;
sin_freq = freq_L;
while true
    w = 2*pi*sin_freq;
    A = 1;
    dot1 = linspace(1, L*2, L*2);
    dot2 = zeros([1, L*2]);
    t_dot = [dot2 ; dot1];

%     figure();
%     plot(dot1, sin(w*dot1*Ts));
%     title("input");

    y = zeros([q, L]);
    for k=1:L
        step = 1;
        for tau=1:L
           u = sin(w*t_dot*Ts);
           y(:, k) = y(:, k) + Yi(:, step:step+m-1) * u(:, k + tau); 
           step = step + m;
        end
    end
    sin_freq = sin_freq * time_value;
    
    [Am1, Ph1] = getSineAmplitudeAndPhase(y(1, :));
    [Am2, Ph2] = getSineAmplitudeAndPhase(y(2, :));
    collect_Am1(1, cnt) = 20*log10(Am1);
    collect_Ph1(1, cnt) = Ph1;
    collect_Am2(1, cnt) = 20*log10(Am2);
    collect_Ph2(1, cnt) = Ph2;
    cnt = cnt + 1;
    if cnt > freq_size
        break
    end
end


freq_x = linspace(freq_L, freq_H, freq_size);

figure();
subplot(2,1,1)
plot(freq_x, collect_Am1);
title("AM1 by input2");
xlabel("Frequency (Hz)");
ylabel("db");
grid on;
set(gca, 'XScale', 'log');
subplot(2,1,2)
plot(freq_x, collect_Ph1);
title("Ph1 by input2");
xlabel("Frequency (Hz)");
ylabel("度");
grid on;
set(gca, 'XScale', 'log');

figure();
subplot(2,1,1)
plot(freq_x, collect_Am2);
title("AM2 by input2");
xlabel("Frequency (Hz)");
ylabel("db");
grid on;
set(gca, 'XScale', 'log');
subplot(2,1,2)
plot(freq_x, collect_Ph2);
title("Ph2 by input2");
xlabel("Frequency (Hz)");
ylabel("度");
grid on;
set(gca, 'XScale', 'log');


