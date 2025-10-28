%Calculate 2d TM RCS 
function [RCS] = RCS_2d(Ez,lambda,Totalnx,Totalny,dx,dy)
c=3e8;             % Speed of light
k0=2*pi/lambda;         % Wavenumber
omega0=k0*c;       % Angular frequency
mu0 = 4*pi*10^(-7); % Vacuum permeability [H/m]
epsilon0=1/(4*pi*9*10^9);
Z0=sqrt(mu0/epsilon0);
%% Calculate the magnetic field using Ez
% Yee grid Hx (i, j) in the (i, j + 0.5), Hy (i, j) in the (i + 0.5, j)
Hx = zeros(Totalnx, Totalny-1);
Hy = zeros(Totalnx-1, Totalny);
% Calculate staggered Hx and Hy
% Calculate Hx
Hx(:, 1:Totalny-1) = 1 / (1i * omega0 * mu0) * (Ez(:, 2:Totalny) - Ez(:, 1:Totalny-1)) / dy;

% Calculate Hy
Hy(1:Totalnx-1, :) = -1 / (1i * omega0 * mu0) * (Ez(2:Totalnx, :) - Ez(1:Totalnx-1, :)) / dx;
% Interpolate to the Ez grid points
Hx_c = zeros(Totalnx, Totalny);
Hy_c = zeros(Totalnx, Totalny);
Hx_c(:,2:Totalny-1) = 0.5 * (Hx(:,1:end-1) + Hx(:,2:end));     % Hx interpolation
Hy_c(2:Totalnx-1,:) = 0.5 * (Hy(1:end-1,:) + Hy(2:end,:));     % Hy interpolation
%% The equivalent surface
NPML=8;
buffer = NPML + 10; % Keep 10 grids away from the PML
boundary_pts = [];
normals = [];
% Bottom boundary
for i = buffer:Totalnx-buffer
    boundary_pts(end+1,:) = [i, buffer];
    normals(end+1,:) = [0, -1];
end
% Top boundary
for i = buffer:Totalnx-buffer
    boundary_pts(end+1,:) = [i, Totalny-buffer];
    normals(end+1,:) = [0, 1];
end
% Left boundary
for j = buffer+1:Totalny-buffer-1
    boundary_pts(end+1,:) = [buffer, j];
    normals(end+1,:) = [-1, 0];
end
% Right boundary
for j = buffer+1:Totalny-buffer-1
    boundary_pts(end+1,:) = [Totalnx-buffer, j];
    normals(end+1,:) = [1, 0];
end

X = (0:Totalnx-1)*dx;
Y = (0:Totalny-1)*dy;
%% Calculate the equivalent electromagnetic current
N = size(boundary_pts,1);
J = zeros(N,1);
Mx = zeros(N,1);
My = zeros(N,1);

for k = 1:N
    i = boundary_pts(k,1);
    j = boundary_pts(k,2);
    n = normals(k,:);
    
    Ez_val = Ez(i,j);

    Hx_val = Hx_c(i,j);
    Hy_val = Hy_c(i,j);

    J(k)  = n(2)*Hx_val - n(1)*Hy_val;
    Mx(k) = -n(2)*Ez_val;
    My(k) =  n(1)*Ez_val;
end
%% Far-field radiation calculation
theta = linspace(0, 2*pi, 361);
E_far = zeros(size(theta));
for m = 1:length(theta)
    phi = theta(m);
    Eaccum = 0;
    for k = 1:N
        i = boundary_pts(k,1);
        j = boundary_pts(k,2);
        x = X(i);
        y = Y(j);
        
        phase = exp(1i * k0 * (x * cos(phi) + y * sin(phi)));
        Eaccum = Eaccum + ( -Z0 * J(k) - Mx(k)*sin(phi) + My(k)*cos(phi) ) * phase;
    end
    E_far(m) = Eaccum * dx;  
end
rcs=(k0/4)*abs(E_far).^2/lambda;
RCS = 10*log10(rcs); 
end