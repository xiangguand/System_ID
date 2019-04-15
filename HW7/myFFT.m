function [fft_signal] = myFFT(signal)
dim_signal = size(signal);
dim_row = dim_signal(1);
dim_col = dim_signal(2);
N = dim_col;
% Do FFT
fft_signal = zeros([dim_row, dim_col]);

for k=1:N
    temp = zeros([dim_row, 1]);
    for n = 1:N
        temp = temp + signal(n) * exp(-i*2*pi*k*n/N);
    end
    fft_signal(:, k) = temp;
end
fft_signal = abs(fft_signal);

end

