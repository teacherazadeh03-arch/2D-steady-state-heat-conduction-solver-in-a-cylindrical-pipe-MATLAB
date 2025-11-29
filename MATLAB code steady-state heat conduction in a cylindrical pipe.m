clc; clear;

%% ---- SOLVER ----
L    = 3.0;      ri = 0.05;  ro = 0.10;
Nr   = 80;       Nz = 200;
r    = linspace(ri, ro, Nr).';   dr = r(2)-r(1);
z    = linspace(0, L, Nz);       dz = z(2)-z(1);

k    = 16;       % W/m-K
Thot = 120;      Tcold = 30;

% Outer convection (Robin)
h    = 15;       Tinf = 25;

% Inner wall insulated
use_inner_convection = false; h_in = 200; Tin = 80;

% Initial guess + axial Dirichlet BCs
T = repmat(linspace(Thot, Tcold, Nz), Nr, 1);
T(:,1) = Thot; T(:,end) = Tcold;

% Cylindrical Laplace coefficients
Arp = (1/dr^2) + 1./(2*r*dr);
Arm = (1/dr^2) - 1./(2*r*dr);
Az  = 1/dz^2;
Ap  = (2/dr^2) + (2/dz^2);

max_iter = 20000; tol = 1e-6;
for iter = 1:max_iter
    T_old = T;

    % Interior
    for i = 2:Nr-1
        for j = 2:Nz-1
            T(i,j) = ( Arp(i)*T(i+1,j) + Arm(i)*T(i-1,j) + Az*(T(i,j+1) + T(i,j-1)) ) / Ap;
        end
    end

    % Inner wall
    if use_inner_convection
        for j = 2:Nz-1
            T(1,j) = ((k/dr)*T(2,j) + h_in*Tin) / (k/dr + h_in);
        end
    else
        T(1,2:Nz-1) = T(2,2:Nz-1); % insulated
    end

    % Outer wall (Robin)
    for j = 2:Nz-1
        T(end,j) = ((k/dr)*T(end-1,j) + h*Tinf) / (k/dr + h);
    end

    % Axial Dirichlet
    T(:,1) = Thot; T(:,end) = Tcold;

    if max(abs(T(:)-T_old(:))) < tol
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end
end

%% ---- 3D CYLINDRICAL VISUALIZATIONS ----

% Data on the outer wall
T_wall = T(end, :);                             % T at r = ro
% Axial gradient on wall
dTdz_wall = zeros(1, Nz);
dTdz_wall(2:Nz-1) = (T_wall(3:Nz) - T_wall(1:Nz-2)) / (2*dz);
dTdz_wall(1)      = (T_wall(2)   - T_wall(1))   / dz;
dTdz_wall(end)    = (T_wall(end) - T_wall(end-1)) / dz;

% Radial gradient (for heat flux) on wall
dTdr_wall = (T(end,:) - T(end-1,:)) / dr;       % ∂T/∂r at r = ro
q_wall    = -k * dTdr_wall;

% Build a cylinder mesh at r = ro
nTheta = 120;                                   % angular resolution for smooth look
theta  = linspace(0, 2*pi, nTheta);
[Theta, Zc] = meshgrid(theta, z);               % (Nz x nTheta)

X = ro * cos(Theta);
Y = ro * sin(Theta);
Z = Zc;

% Map 1D wall arrays (size 1xNz) to (Nz x nTheta) for color data
Twall_map    = repmat(T_wall(:),    1, nTheta);
dTdz_map     = repmat(dTdz_wall(:), 1, nTheta);
qwall_map    = repmat(q_wall(:),    1, nTheta);

% --- 3D Cylinder colored by Temperature on the wall ---
figure;
surf(X, Y, Z, Twall_map, 'EdgeColor','none'); axis equal;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Outer Wall Temperature on 3D Pipe Surface');
colormap(turbo); colorbar; view(135,25);
camlight headlight; lighting gouraud;

% --- 3D Cylinder colored by Axial Gradient (∂T/∂z) on the wall ---
figure;
surf(X, Y, Z, dTdz_map, 'EdgeColor','none'); axis equal;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Axial Thermal Gradient (\partialT/\partialz) on Outer Wall');
colormap(parula); colorbar; view(-35,25);
camlight headlight; lighting gouraud;

% --- 3D Cylinder colored by Heat Flux (q_r) on the wall (optional) ---
figure;
surf(X, Y, Z, qwall_map, 'EdgeColor','none'); axis equal;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title('Radial Heat Flux q_r = -k \partialT/\partialr on Outer Wall');
colormap(hot); colorbar; view(35,25);
camlight headlight; lighting gouraud;

%% (Optional) Diagnostic r–z surface for the full field
[Rg, Zg] = meshgrid(r, z);
figure;
surf(Zg, Rg, T.', 'EdgeColor','none'); view(135,30);
xlabel('z (m)'); ylabel('r (m)'); zlabel('T (°C)');
title('Temperature Field in Pipe Wall (r–z)');
colormap(turbo); colorbar; shading interp; grid on;
