function [A, P_now, K, sigma, e_hat] = ExpWinRLS(x, y, P_before, A, lambda)
    [r, c] = size(x);
    sigma = lambda + (x * P_before * x');  % 1x1
    K = P_before * x' * (sigma^-1);
    P_now = (eye(r) - K * x) * (lambda^-1) * P_before;
    e_hat = y - x*A;      % true output - estimate output = estimation of error
    A = A + K * e_hat;    % update the A parameters
