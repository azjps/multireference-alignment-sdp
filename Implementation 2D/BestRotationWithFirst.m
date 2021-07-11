% returns best l_2 shift with first image

function avg_img = BestRotationWithFirst(Y, sigma, DEBUG)
    [R,L,N] = size(Y);

    shifts = zeros(N,1);
    avg_img = Y(:,:,1);
    for n = 2:N
        [~, shifts(n)] = ShiftLessDist(Y(:,:,1), Y(:,:,n));
        avg_img = avg_img + RotateImage(Y(:,:,n), shifts(n), 0);
    end

    avg_img = avg_img / N;
end