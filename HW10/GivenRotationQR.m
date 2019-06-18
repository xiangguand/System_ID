function [Q, R] = GivenRotationQR(A)
    [m, n] = size(A);
    Q = zeros([m, m]);
    R = A;
    for ni=1:n-1+(m-n)
        for mi=ni+1:m
            % eliminate A(mi, ni)
            if ni < n
                c = R(ni, ni);
                s = R(mi, ni);
            else
                c = R(ni, n);
                s = R(mi, n);
            end
            r = sqrt(c^2 + s^2);
            if r ~= 0
                c = c/r;
                s = s/r;
            else
                c = 1;
                s = 1;
            end
            q = eye(m);
            q(ni, ni) = c;
            q(mi, mi) = c;
            q(ni, mi) = s;
            q(mi, ni) = -s;
            R = q*R;
            if sum(sum(Q)) == 0
                Q = q';
            else
                Q = Q*q';
            end
        end
    end
end

