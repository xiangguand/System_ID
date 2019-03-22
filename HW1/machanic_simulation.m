clc, clear, close all

config_file = fopen('example_machanic_n2.json');
save_text = ''
while 1
    tline = fgetl(config_file);
    if ~ischar(tline), break, end
    disp(tline);
    if tline ~= ""
        if tline(1) ~= '#'
            save_text = strcat(save_text, tline);
        end
    end
end
fclose(config_file);
data = jsondecode(save_text);

n = data.test_sample.number;
unit = data.test_sample.unit;
input_vect = data.test_sample.input;
mass = data.test_sample.mass;
damping = data.test_sample.damping;
spring = data.test_sample.spring;
mass_sz = size(mass);
damping_sz = size(damping);
spring_sz = size(spring);
para_sz_vect = [mass_sz(1) damping_sz(1) spring_sz(1)];
if ~(sum(para_sz_vect==n)==3)
    disp('Config file parameter numbers error !!!');
    disp('Enter Ctrl+C to break while loop !!!');
    while(1)
    end
else
    disp('Load config file parameters okay');
end

syms s;
para_matrix = sym('s', [n, n]);   % using symbolic tool
for i=1:n
    for j=1:n
        para_matrix(i,j) = 0;      % initialize
    end
end

if n==1
    para_matrix(1) = mass(1)*s^2 + damping(1)*s + spring(1);
elseif n==2
    para_matrix(1,:) = [mass(1)*s^2+(damping(1)+damping(2))*s+spring(1)+spring(2), -(damping(2)*s+spring(2))];
    para_matrix(2,:) = [-(damping(2)*s+spring(2)), mass(2)*s^2+damping(2)*s+spring(2)];
else
    para_matrix(1, 1:2) = [mass(1)*s^2+(damping(1)+damping(2))*s+spring(1)+spring(2), -(damping(2)*s+spring(2))]; % first one
    for k=2:n-1
        para_matrix(k, k-1:k+1) = [-(damping(k)*s+spring(k)), mass(k)*s^2+(damping(k)+damping(k+1))*s+spring(k)+spring(k+1), -(damping(k+1)*s+spring(k+1))];
    end
    para_matrix(n, n-1:n) = [-(damping(n)*s+spring(n)), mass(n)*s^2+damping(n)*s+spring(n)]; % last one
end

para_matrix

% using Cramer's Rule
delta = det(para_matrix);
temp_para = para_matrix;
X_n = sym('s', [1, n]);
G_n = sym('s', [1, n]);
for k=1:n
    temp_para(:, k) = input_vect;
    X_n(k) =  simplify(det(temp_para)/delta);
    temp_para = para_matrix;
    display(strcat('G_', num2str(k)));
    G_n(k) = X_n(k) / input_vect(k)
    
    % T: transfer function
    % A¡BB¡BC¡BD are the state space parameters
    [T, A, B, C, D] = symArray2StateSpace(G_n(k))

    sampling_time = 0.001;
    T_discrete = c2d(T, sampling_time, 'zoh')
    figure();
    step(T);
    title('Continuous time step');

    figure();
    step(T_discrete)
    title('Discrete time step (ZOH)');
end







