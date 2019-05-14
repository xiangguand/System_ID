function [B] = CGS(X, normalize_en)
    [m, n] = size(X);
    % m: vector dimension
    % n: how many the vectors have.
    B = X;
    for j=1:n
        for i=1:j-1
            B(:, j) = B(:, j) - dot(B(:, j), B(:, i)) / dot(B(:, i), B(:, i)) * B(:, i);
        end
        % Do nomailize
        if normalize_en
            B(:, j) = B(:, j) / norm(B(:, j));
        end
    end
end

