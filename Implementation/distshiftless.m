% Returns the sum of the minimal distance between a rotation of x and z.
function sum_dist = distshiftless(x,Z)

    [L,N] = size(Z);

    sum_dist = 0;
    for n = 1:N
        z = Z(:,n);
        dist = norm(x-z);

        for k=1:L-1 
              z = circshift(z,1);
              if norm(x-z)<dist
                   dist = norm(x-z); 
              end
        end
        sum_dist = sum_dist + dist;
    end
end