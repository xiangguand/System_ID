function [signal] = myIFFT(fft_signal)
dim_signal = size(fft_signal);
dim_row = dim_signal(1);
dim_col = dim_signal(2);
N = dim_col;
% Do FFT
signal = zeros([dim_row, dim_col]);

for k=1:N
    temp = zeros([dim_row, 1]);
    for n = 1:N
        temp = temp + fft_signal(n) * exp(i*2*pi*k*n/N);
    end
    signal(:, k) = temp/N;
end
signal = abs(signal);

end