function [pseudoInv] = pseudo_inverse(mat)
    if det(mat) < 0.0000000001
        pseudoInv = 0;
    else
        pseudoInv = mat^-1;
    end
        
end

