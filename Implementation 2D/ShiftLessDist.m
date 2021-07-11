% Returns the sum of the minimal distance between a rotation of Y and Z.
function [sum_dist, min_k] = ShiftLessDist(Y,Z)

    [~, L] = size(Z);
    
    dist = @(k) norm(double(Y) - RotateImage(Z, k, 0), 'fro');
       
    [sum_dist, min_k] = min(arrayfun(dist, 1:L));
    
end