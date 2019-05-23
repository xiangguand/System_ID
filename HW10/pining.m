function [pining_vector] = pining(row_size, column_size)
    pining_vector = zeros([row_size, column_size]);
    pining_vector(1, :) = 1;
end

