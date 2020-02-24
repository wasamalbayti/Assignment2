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
            F(n) = 1;
            
        elseif j == 1
            G(n,:) = 0;
            G(n,n) = 1;
            
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
            
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

Emap = zeros(ny,nx,1);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny; 
        Emap(j,i) = V(n);
    end
end

figure(1)
surf(Emap)
title("Numerical Method Voltage Surface Plot")
xlabel("X")
ylabel("Y")
zlabel("Voltage")

% Now to plot the analytical solution
a = ny;
b = nx/2;
[x,y] = meshgrid((linspace(-nx/2,nx/2,20)),(linspace(0,ny,ny)));
Emap2 = sparse(ny,nx);

% now to create a simulation for anayltical solution
for n = 1:2:500
    Emap2 = Emap2+(4./pi).*((cosh(n*pi*x/a).*sin(n*pi*y/a))./(n*cosh(n*pi*b/a)));
    figure(2)
    surf(Emap2)
    title("Analytical Method Voltage Surface Plot")
    xlabel("X")
    ylabel("Y")
    zlabel("Voltage")
end
