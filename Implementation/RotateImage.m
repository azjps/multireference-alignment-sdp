
function Z = RotateImage(Y, L, sigma)

    [R C] = size(Y);
    Z = zeros(size(Y));
    for i = 1:R
        Z(i,:) = circshift(Y(i,:)', L)'; 
        if sigma ~= 0
            for j = 1:C
                Z(i,j) = Z(i,j) + round(normrnd(0, sigma)); 
            end
        end
    end

end