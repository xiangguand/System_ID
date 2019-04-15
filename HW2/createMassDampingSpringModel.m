function [A, B, C, D, para_struct] = createMassDampingSpringModel(json_config_file)

config_file = fopen(json_config_file);
save_text = '';
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
    error('Config file parameter numbers error !!!');
else
    disp('Load config file parameters okay');
end

%%% Save parameters to structure type
para_struct.n = n;
para_struct.input_sz = n;
para_struct.output_sz = n;
para_struct.state_sz = n*2;
para_struct.unit = unit;
para_struct.input_vect = input_vect;
para_struct.mass = mass;
para_struct.damping = damping;
para_struct.spring = spring;

syms s;
para_matrix = sym('s', [n, n+1]);   % using symbolic tool
para_matrix(:, end) = ones([1, n]); % input parameters
for i=1:n
    for j=1:n
        para_matrix(i,j) = 0;      % initialize
    end
end

% initialize state space parameters
A = zeros([n*2, n*2]);
B = zeros([n*2, n]);
C = zeros([n, n*2]);
D = zeros([n, n]);

if n==1
    para_matrix(1) = mass(1)*s^2 + damping(1)*s + spring(1);
elseif n==2
    para_matrix(1,1:end-1) = [mass(1)*s^2+(damping(1)+damping(2))*s+spring(1)+spring(2), -(damping(2)*s+spring(2))];
    para_matrix(2,1:end-1) = [-(damping(2)*s+spring(2)), mass(2)*s^2+damping(2)*s+spring(2)];
else
    para_matrix(1, 1:2) = [mass(1)*s^2+(damping(1)+damping(2))*s+spring(1)+spring(2), -(damping(2)*s+spring(2))]; % first one
    for k=2:n-1
        para_matrix(k, k-1:k+1) = [-(damping(k)*s+spring(k)), mass(k)*s^2+(damping(k)+damping(k+1))*s+spring(k)+spring(k+1), -(damping(k+1)*s+spring(k+1))];
    end
    para_matrix(n, n-1:n) = [-(damping(n)*s+spring(n)), mass(n)*s^2+damping(n)*s+spring(n)]; % last one
end
para_matrix = para_matrix./mass     % let highest s be 1

%%% Save para_matrix to structure type
para_struct.para_matrix = para_matrix;

% create A, B, C, D state space
for i=1:n*2
    if mod(i, 2) ~= 0
       A(i, i+1) = 1;
       %%% C matrix
       C((i-1)/2+1, i) = 1;
    else
        %%% A matrix
        temp_arr = zeros([1, n*2]);
        if i <= 4
            pos = 1;
        else
            pos = pos + 2;
        end
        p = pos;
        for j=1:n
            if para_matrix(i/2, j) ~= 0
               t = sym2poly(para_matrix(i/2, j));
               t = t(:, end-1:end);
               temp_arr(1, p:p+1) = -t(:, end:-1:1);
               p = p + 2;
            end
        end
        A(i, :) = temp_arr;
        %%% B matrix
        B(i, i/2) = para_matrix(i/2, end);
    end
    
end


end

