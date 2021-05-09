% Constuction of a quadrilateral 3D mesh

% Domain and element size description
Lx = 3;Ax = 1/2;
Ly = 1;Ay = 1/2;
Lz = 1;Az = 1/2;

mx = Lx/Ax;
my = Ly/Ay;
mz = Lz/Az;

X = zeros((mx+1)*(my+1)*(mz+1), 3);
R = zeros(mx*my*mz,8);

% Construction
cont = 1;
for k = 1:(mz+1)
    for j = 1:(my+1)
        for i = 1:(mx+1)
            X(cont,:) = [Ax*(i-1), Ay*(j-1), Az*(k-1)];
            cont = cont+1;
        end
    end
end

cont = 1;
for m = 0:mz-1
    for l = 0:my-1
        for k = 1:mx
            i = k + l*(mx+1) + m*(mx+1)*(my+1);
            R(cont,:) = [i, i+1, mx+2+i, mx+1+i, (mx+1)*(my+1)+i, (mx+1)*(my+1)+i+1, (mx+1)*(my+2)+i+1, (mx+1)*(my+2)+i];
            cont = cont+1;
        end
    end
end

save('mesh9.mat', 'X', 'R', 'Lx', 'Ly', 'Lz', 'Ax', 'Ay', 'Az')



