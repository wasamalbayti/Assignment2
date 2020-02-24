% this section will involve a bottleneck and variosu plots are required
% Dimensions of matrices, using ratio of 3/2 for L/W
nx = 20;
ny = 30; 

% G matrix so the solution will be similifed to a matrix form i.e Ax = B
G = sparse(nx*ny);
F = zeros(1,nx*ny);

% Set up the sigma matrix 
sigma = zeros(nx,ny); % initalize the sigma matrix to zeros at first
for i = 1:nx
    for currdensity = 1:ny
        if i > (2/5)*nx && i < (3/5)*nx && (currdensity < (2/5)*ny||currdensity > (3/5)*ny)
            sigma(i, currdensity) = 10^-2;
        else
            sigma(i, currdensity) = 1;
        end
    end
end

for x = 1:nx
    for y = 1:ny
        n = y + (x-1)*ny;
        nxp = y + (x+1-1)*ny;
        nxx = y + (x-1-1)*ny;
        nyp = y + 1 + (x-1)*ny;
        nnegy = y - 1 + (x-1)*ny;
        if x == 1       
            G(n, :) = 0;
            G(n, n) = 1;
            F(n) = 1;
        elseif x == nx
            G(n, :) = 0;
            G(n, n) = 1;
            F(n) = 0;
        elseif y == 1
            G(n, nxp) = (sigma(x+1, y) + sigma(x,y))*(1/2);
            G(n, nxx) = (sigma(x-1, y) + sigma(x,y))*(1/2);
            G(n, nyp) = (sigma(x, y+1) + sigma(x,y))*(1/2);            
            G(n, n) = -(G(n,nxp)+G(n,nxx)+G(n,nyp));
        elseif y == ny
            G(n, nxp) = (sigma(x+1, y) + sigma(x,y))*(1/2);
            G(n, nxx) = (sigma(x-1, y) + sigma(x,y))*(1/2);
            G(n, nnegy) = (sigma(x, y-1) + sigma(x,y))*(1/2);
            G(n, n) = -(G(n,nxp)+G(n,nxx)+G(n,nnegy));    
        else
            G(n, nyp) = (sigma(x, y+1) + sigma(x,y))*(1/2);
            G(n, nnegy) = (sigma(x, y-1) + sigma(x,y))*(1/2);
            G(n, nxp) = (sigma(x+1, y) + sigma(x,y))*(1/2);
            G(n, nxx) = (sigma(x-1, y) + sigma(x,y))*(1/2);
            G(n, n) = -(G(n,nxp)+G(n,nxx)+G(n,nyp)+G(n,nnegy));
        end
    end
end
Voltage = G\F';

%flip the axis 
Emap = zeros(ny, nx, 1);
for i = 1:nx
    for currdensity = 1:ny
        n = currdensity + (i-1)*ny;
        Emap(currdensity,i) = Voltage(n);
    end
end

% electric field
[Ex, Ey] = gradient(Emap);

% curent desnity calculation
currdensityx = sigma'.*Ex;
currdensityy = sigma'.*Ey;
currdensity = sqrt(currdensityx.^2 + currdensityy.^2);

figure(1)
title("Sigma surface plot")
surf(sigma);
xlabel("X")
ylabel("Y")
zlabel("Sigma")

figure(2)
surf(Emap)
title("Voltage surface plot")
xlabel("X")
ylabel("Y")
zlabel("V")

figure(3)
surf(-Ex)
title("X surface plot")
xlabel("X")
ylabel("Y")
zlabel("V")

figure(4)
surf(-Ey)
title("y surface plot")
xlabel("X")
ylabel("Y")
zlabel("V")


figure(5)
surf(currdensity)
title("Curent Density Plot")
xlabel("X")
ylabel("Y")
zlabel("Current Density")


