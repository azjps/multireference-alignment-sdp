function understandingDual

%L =2;
% N =3;

goon = 1;

L = zeros(3,3);

while goon

L(1,2) = random('normal',0,1);
L(1,2) = -abs(L(1,2));
L(2,1) = L(1,2);
L(1,3) = random('normal',0,1);
L(1,3) = -abs(L(1,3));
L(3,1) = L(1,3);
L(1,1) = abs(L(1,2)) + abs(L(1,3));
L(2,3) = random('normal',0,1);
L(2,3) = -abs(L(2,3));
L(3,2) = L(2,3);
L(2,2) = abs(L(2,1)) + abs(L(2,3));
L(3,3) = abs(L(3,1)) + abs(L(3,2));

La = L

L(1,2) = random('normal',0,1);
L(1,2) = -abs(L(1,2));
L(2,1) = L(1,2);
L(1,3) = random('normal',0,1);
L(1,3) = -abs(L(1,3));
L(3,1) = L(1,3);
L(1,1) = abs(L(1,2)) + abs(L(1,3));
L(2,3) = random('normal',0,1);
L(2,3) = -abs(L(2,3));
L(3,2) = L(2,3);
L(2,2) = abs(L(2,1)) + abs(L(2,3));
L(3,3) = abs(L(3,1)) + abs(L(3,2));

Lb = L

if La(1,1)>Lb(1,1)
    if La(2,2)>Lb(2,2)
        if La(3,3)>Lb(3,3)
            goon = 0;
        end
    end
end

end

I2 = eye(2);
O2 = [0,1;1,0];

ALa = kron(La,I2);
ALb = kron(Lb,O2);
Q = ALa + ALb

eigALa = eig(ALa)'
eigALb = eig(ALb)'
eigQ = eig(Q)'

end