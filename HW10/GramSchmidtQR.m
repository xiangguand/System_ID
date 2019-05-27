function [Q, R] = GramSchmidtQR(A)
    [m, n] = size(A);
    Q = zeros([m, n]);
    R = zeros([n, n]);
    q_hat = 1;
    for ni=1:n
       if ni==1
           q = A(:, ni)/norm(A(:, ni));
       else
          for i=1:ni
             if i==1
                q_hat = eye(m) - Q(:, i)*Q(:, i)';
             else
                 q_hat = q_hat * (eye(m)-Q(:, i)*Q(:, i)');
             end
          end
          q_hat = q_hat * A(:, ni);
          q = q_hat/norm(q_hat);
       end
       % update R
       for j=1:ni
          if j==ni
              if ni==1
                 R(j, ni) = norm(A(:, ni));
              else
                 R(j, ni) = norm(q_hat);   
              end
          else
              R(j, ni) = Q(:, j)'*A(:, ni);
          end
       end
       Q(:, ni) = q;
    end
end

