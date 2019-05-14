function [V] = MGS(X)
    [m, n] = size(X);
    % m: vector dimension
    % n: how many the vectors have.
    V = X;
    for j=1:n
        % Do normalize
        temp_v = V(:, j) / norm(V(:, j)');
        for k=j+1:n   % calculate other vectors
            V(:, k) = V(:, k) - temp_v'*V(:, k)*temp_v;
        end
    end
end

