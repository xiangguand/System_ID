function [projection_vector] = projection(vector, R)
    % vector_i : Length X dimension
    % project vector to R
    projection_vector = R * pseudo_inverse(R' * R) * R' * vector;
end

