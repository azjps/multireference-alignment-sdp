function optshift = optshiftfrom2to1(x,z)

    dist = norm(x-z);
    optshift = 0;
    
    L = length(x);

    for k=1:L-1 
          z = circshift(z,1);
          if norm(x-z)<dist
               dist = norm(x-z);
               optshift = k;
          end
    end
    
end