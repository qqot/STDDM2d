function [Subdomain] = AssemblyofEqu(Subdomain,k0,dx,dy,Ndomain)
%Assemble system matrix
c=3e8;             % Speed of light
eps0=8.854e-12;
omega0=k0*c;       % Angular frequency
PMLnx=8;
PMLny=8;
Totalnx=Subdomain.Totalsubnx;
Totalny=Subdomain.Totalsubny;

Subdomain.Eb=zeros(Totalnx*Totalny,1); %backgroud electric field
%% PML parameters
sigma_power=3.7;
R0=100;
eta=377;
sigmax_max=(sigma_power+1)*R0/eta/2/(PMLnx*dx)/k0/(sqrt(Ndomain));
sigmax_max=0.3;

% Initialize all sigma arrays
sigmax_int = zeros(1, Totalnx);
sigmax_hfint = zeros(1, Totalnx+1);
sigmay_int = zeros(1, Totalny);
sigmay_hfint = zeros(1, Totalny+1);

% Left PML layers
pml_indices = 1:PMLnx;
sigmax_int(pml_indices) = sigmax_max/dx * ((2*(PMLnx:-1:1)-1)/(2*PMLnx)).^sigma_power;
sigmax_hfint(1:PMLnx+1) = sigmax_max/dx * ((PMLnx:-1:0)/PMLnx).^sigma_power;

% Right PML layers
pml_indices = Totalnx-PMLnx+1:Totalnx;
sigmax_int(pml_indices) =  sigmax_max/dx * ((2*(1:PMLnx)-1)/(2*PMLnx)).^sigma_power;
sigmax_hfint(Totalnx-PMLnx+1:Totalnx+1) =  sigmax_max/dx * ((0:PMLnx)/PMLnx).^sigma_power;

% Bottom PML layers
pml_indices = 1:PMLny;
sigmay_int(pml_indices) =  sigmax_max/dy * ((2*(PMLny:-1:1)-1)/(2*PMLny)).^sigma_power;
sigmay_hfint(1:PMLny+1) = sigmax_max/dy * ((PMLny:-1:0)/PMLny).^sigma_power;

% Top PML layers
pml_indices = Totalny-PMLny+1:Totalny;
sigmay_int(pml_indices) =  sigmax_max/dy * ((2*(1:PMLny)-1)/(2*PMLny)).^sigma_power;
sigmay_hfint(Totalny-PMLny+1:Totalny+1) =  sigmax_max/dy * ((0:PMLny)/PMLny).^sigma_power;

kx_int=1-1i*sigmax_int/omega0/eps0;
ky_int=1-1i*sigmay_int/omega0/eps0;
kx_hfint=1-1i*sigmax_hfint/omega0/eps0;
ky_hfint=1-1i*sigmay_hfint/omega0/eps0;
%% FDFD Matrix Formulation
epsilon=Subdomain.eps;

N = Totalnx * Totalny;

% Preallocate arrays with maximum needed size (5 values per grid point except edges)
max_elements = 5*N;
m = zeros(1, max_elements);
n = zeros(1, max_elements);
s = zeros(1, max_elements);
% Create source vector
src = zeros(N, 1);
% Counter to keep track of the actual number of elements added
count = 0;

% Calculate row and column indices for each grid point
[j_grid, i_grid] = meshgrid(1:Totalny, 1:Totalnx);
i_grid = i_grid(:);
j_grid = j_grid(:);
indices = (j_grid-1)*Totalnx + i_grid;

% Interior points (efficient bulk processing)
interior_mask = (i_grid > 1) & (i_grid < Totalnx) & (j_grid > 1) & (j_grid < Totalny);
interior_indices = indices(interior_mask);
i_int = i_grid(interior_mask);
j_int = j_grid(interior_mask);

% Process all interior points using true vectorization
num_interior = length(interior_indices);
interior_eps = epsilon(sub2ind(size(epsilon), i_int, j_int));

% Calculate positions for incident field
xx_interior = (i_int-1)*dx+Subdomain.subx(1);
yy_interior = (j_int-1)*dy+Subdomain.suby(1);

% Diagonal elements (center)
diag_start = count + 1;
m(diag_start:diag_start+num_interior-1) = interior_indices;
n(diag_start:diag_start+num_interior-1) = interior_indices;
s(diag_start:diag_start+num_interior-1) = k0^2*interior_eps.' - ...
    1./dx^2./kx_int(i_int).*(1./kx_hfint(i_int)+1./kx_hfint(i_int+1)) - ...
    1./dy^2./ky_int(j_int).*(1./ky_hfint(j_int)+1./ky_hfint(j_int+1));

% Left neighbor elements
left_start = diag_start + num_interior;
m(left_start:left_start+num_interior-1) = interior_indices;
n(left_start:left_start+num_interior-1) = interior_indices - 1;
s(left_start:left_start+num_interior-1) = 1./dx^2./kx_int(i_int)./kx_hfint(i_int);

% Right neighbor elements
right_start = left_start + num_interior;
m(right_start:right_start+num_interior-1) = interior_indices;
n(right_start:right_start+num_interior-1) = interior_indices + 1;
s(right_start:right_start+num_interior-1) = 1./dx^2./kx_int(i_int)./kx_hfint(i_int+1);

% Bottom neighbor elements
bottom_start = right_start + num_interior;
m(bottom_start:bottom_start+num_interior-1) = interior_indices;
n(bottom_start:bottom_start+num_interior-1) = interior_indices - Totalnx;
s(bottom_start:bottom_start+num_interior-1) = 1./dy^2./ky_int(j_int)./ky_hfint(j_int);

% Top neighbor elements
top_start = bottom_start + num_interior;
m(top_start:top_start+num_interior-1) = interior_indices;
n(top_start:top_start+num_interior-1) = interior_indices + Totalnx;
s(top_start:top_start+num_interior-1) = 1./dy^2./ky_int(j_int)./ky_hfint(j_int+1);

count = top_start + num_interior - 1;


% Process boundary points (PML edges) using vectorization
% Left and Right boundaries
left_i = ones(1,Totalny);
right_i = Totalnx*ones(1,Totalny);
boundary_j = 1:Totalny;

% Combine left and right boundary indices
boundary_i = [left_i, right_i];
boundary_j = [boundary_j, boundary_j];
boundary_indices = boundary_i + (boundary_j-1)*Totalnx;

% Vectorized calculation for left and right PML edges
num_lr_boundary = length(boundary_indices);
boundary_start = count + 1;
m(boundary_start:boundary_start+num_lr_boundary-1) = boundary_indices;
n(boundary_start:boundary_start+num_lr_boundary-1) = boundary_indices;
s(boundary_start:boundary_start+num_lr_boundary-1) = k0^2 - ...
    1./dx^2./kx_int(boundary_i).*(1./kx_hfint(boundary_i)+1./kx_hfint(boundary_i+1)) - ...
    1./dy^2./ky_int(boundary_j).*(1./ky_hfint(boundary_j)+1./ky_hfint(boundary_j+1));

count = boundary_start + num_lr_boundary - 1;

% Top and Bottom boundaries (excluding corners already counted)
top_bottom_i = 2:Totalnx-1;
bottom_j = ones(1,Totalnx-2);
top_j = Totalny*ones(1,Totalnx-2);

% Combine top and bottom indices
boundary_i = [top_bottom_i, top_bottom_i];
boundary_j = [bottom_j, top_j];
boundary_indices = boundary_i + (boundary_j-1)*Totalnx;

% Vectorized calculation for top and bottom PML edges
num_tb_boundary = length(boundary_indices);
boundary_start = count + 1;
m(boundary_start:boundary_start+num_tb_boundary-1) = boundary_indices;
n(boundary_start:boundary_start+num_tb_boundary-1) = boundary_indices;
s(boundary_start:boundary_start+num_tb_boundary-1) = k0^2 - ...
    1./dx^2./kx_int(boundary_i).*(1./kx_hfint(boundary_i)+1./kx_hfint(boundary_i+1)) - ...
    1./dy^2./ky_int(boundary_j).*(1./ky_hfint(boundary_j)+1./ky_hfint(boundary_j+1));

count = boundary_start + num_tb_boundary - 1;

% Incident field calculation (source term)
incident_mask = (i_int > PMLnx) & (i_int < Totalnx-PMLnx) & ...
    (j_int > PMLny) & (j_int < Totalny-PMLny);
incident_indices = interior_indices(incident_mask);
incident_xx = xx_interior(incident_mask);
incident_eps = interior_eps(incident_mask);
eb0 = Eb_func(incident_xx, k0);
src(incident_indices) = (1-incident_eps).*k0^2.*eb0;
% Trim excess zeros from arrays
m = m(1:count);
n = n(1:count);
s = s(1:count);
A=sparse(m,n,s);
Subdomain.A=A;
Subdomain.b=src;
end