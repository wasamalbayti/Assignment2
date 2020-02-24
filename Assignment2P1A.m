% Dimensions of matrices, using ratio of 3/2 for L/W
nx = 20;
ny = 30; 

% G matrix so the solution will be similifed to a matrix form i.e Ax = B
G = sparse(nx*ny,nx*ny);
F = zeros(nx*ny,1);

% similiar to the PA-5, we will set up the bulk nodes and boundary
% conditions
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 1;
            
        elseif i == nx
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = 0;
            
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = -3;
            G(n,n+1) = 1;
            G(n,n-ny) = 1;
            G(n,n+ny) = 1;
            
        elseif j == ny
            G(n,n) = -3;
            G(n,n-1) = 1;
            G(n,n-ny) = 1;
            G(n,n+ny) = 1;
            
        else 
            G(n,n) = -4;
            G(n,n-1) = 1;
            G(n,n+1) = 1;
            G(n,n-ny) = 1;
            G(n,n+ny) = 1;
        end
    end
end

V = G\F;

Emap = zeros(nx,ny,1);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny; 
        Emap(i,j) = V(n);
    end
end

surf(Emap)
title("Voltage Plot")
xlabel("X")
ylabel("Y")
zlabel("Voltage")

        