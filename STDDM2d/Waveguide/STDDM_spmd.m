%  A source transfer domain decomposition method to calculate
%  topological photonic crystal waveguide
%  The code was written by Zhanwen Wang from Zhejiang University
%  Email: zw_wang@zju.edu.cn

clear;
clc;

%% Parameters
scale=1*1e-3;
% PhC 
R=6*scale;
a0=2.8*R;
Px=112;    %Number of lattice in x-direction
Py=25;     %Number of lattice in y-direction

PMLnx=8;     % Number of PML layers
PMLny=8;     % Number of PML layers
dx=0.2*scale;       % Length of the mesh
dy=0.2*scale;       % Width of the mesh
a = 1663*scale-2*PMLnx*dx;     %Length: (Px-1)*a0-(Py-1)*0.5*a0=1663;             
b=349*scale-2*PMLny*dy;        %Width: (Ny-1)*0.5*sqrt(3)*a0=349
nx = round(a/dx + 1);       % Number of grids (x-direction) in the region of interest
ny =round(b/dy+ 1);         % Number of grids (y-direction) in the region of interest
Totalnx=nx+2*PMLnx;
Totalny=ny+2*PMLny;
Totala=a+2*PMLnx*dx;
Totalb=b+2*PMLny*dy;

c=3e8;             % Speed of light
fre=7.96e9;        %frequency   
lambda0=c/fre;     
k0=2*pi/lambda0;   % Wavenumber
omega0=k0*c;       % Angular frequency
eps0=8.85e-12;  
eps_1=1;
eps_2=11.7;
eps_comsol= Get_eps(a,b,dx,dy,(Py-1)*0.5*a0+R,scale);       %Get eps from comsol
epsilon_center=eps_comsol;
epsilon=ones(Totalnx,Totalny);
i_range = 2:Totalnx-1;
j_range = 2:Totalny-1;

% Average the four surrounding cell centers for each grid point
epsilon(i_range, j_range) = (epsilon_center(i_range-1, j_range-1) + ...
    epsilon_center(i_range-1, j_range) + ...
    epsilon_center(i_range, j_range-1) + ...
    epsilon_center(i_range, j_range)) / 4;

%source
l=dx*4;
x1=(8-1)*0.5*a0+R+a0+24*a0-((Py-1)*0.5*a0+R);
y1=(9-1)*0.5*sqrt(3)*a0-sqrt(3)/2*R;
x2=(9-1)*0.5*a0+R+24*a0-((Py-1)*0.5*a0+R);
y2=(10-1)*0.5*sqrt(3)*a0+sqrt(3)/2*R;
%% Domain Decomposition
%Domain decomposition needs to ensure the topological attributes of each subdomain.
%Domain decomposition can't disrupt trapezoidal shape.
interfacex=[0 10*a0-2*PMLnx*dx  50*a0-2*PMLnx*dx a];   
interfacey=[0 4*0.5*sqrt(3)*a0 20*0.5*sqrt(3)*a0 b];
interfacex(2:end-1)=roundn(interfacex(2:end-1),-3);
interfacey(2:end-1)=roundn(interfacey(2:end-1),-3);
NbrSubdomainx=length(interfacex)-1;NbrSubdomainy=length(interfacey)-1;
Ndomain=NbrSubdomainx*NbrSubdomainy;
Subdomains=cell(Ndomain,1);
idx_matrix = reshape(1:NbrSubdomainx*NbrSubdomainy, NbrSubdomainx, NbrSubdomainy);

% ---------- 启动并行池 ----------
delete(gcp('nocreate'));
parpool(Ndomain);
tic
spmd(Ndomain)
   myID = labindex;
    % Initialize the subdomain structure
    mySub=struct();
    mySub.id=myID;
    % Calculate the position of the subdomain
    [ii, jj] = ind2sub([NbrSubdomainx, NbrSubdomainy], myID);
    mySub.subx = [interfacex(ii), interfacex(ii+1)];
    mySub.suby = [interfacey(jj), interfacey(jj+1)];

    Totalsubnx = round((mySub.subx(2) - mySub.subx(1))/dx + 1 + 2*PMLnx);
    Totalsubny = round((mySub.suby(2) - mySub.suby(1))/dy + 1 + 2*PMLny);
    mySub.Totalsubnx = Totalsubnx;
    mySub.Totalsubny = Totalsubny;

    % Set neighbor information
    neighbors = struct();

    % Bottom neighbor
    if jj > 1
        neighbors.bottom.id = idx_matrix(ii, jj-1);
        neighbors.bottom.send = zeros(1, 9*Totalsubnx);
        neighbors.bottom.recv = zeros(1, 9*Totalsubnx);
    end

    % Top neighbor
    if jj < NbrSubdomainy
        neighbors.top.id = idx_matrix(ii, jj+1);
        neighbors.top.send = zeros(1, 9*Totalsubnx);
        neighbors.top.recv = zeros(1, 9*Totalsubnx);
    end

    % Left neighbor
    if ii > 1
        neighbors.left.id = idx_matrix(ii-1, jj);
        neighbors.left.send = zeros(1, 9*Totalsubny);
        neighbors.left.recv = zeros(1, 9*Totalsubny);
    end

    % Right neighbor
    if ii < NbrSubdomainx
        neighbors.right.id = idx_matrix(ii+1, jj);
        neighbors.right.send = zeros(1, 9*Totalsubny);
        neighbors.right.recv = zeros(1, 9*Totalsubny);
    end

    mySub.neighbors = neighbors;
    fprintf('Worker %d: Complete initialization\n', labindex);
    % Set subdomain parameters and assemble the subsystem matrix
    mySub = SubParameters(mySub, k0, Totalnx, Totalny, epsilon, dx, dy);
    mySub = AssemblyofEqu(mySub, k0, dx, dy,Ndomain);
    mySub.Es=zeros(mySub.Totalsubnx*mySub.Totalsubny,1);
    mySub.Ez=zeros(mySub.Totalsubnx*mySub.Totalsubny,1);
    fprintf('Worker %d: Complete the matrix assembly\n', labindex);
    % Sparse LU factorization
    [mySub.L, mySub.U, mySub.P, mySub.Q] = lu(mySub.A);
    fprintf('Worker %d: Complete factorization\n', labindex);
    % Source
    if x1 > mySub.subx(1) && x1 < mySub.subx(2) && y1 > mySub.suby(1) && y1 < mySub.suby(2)
        x_source = x1-l/2:2*dx:x1+l/2;
        y_source = y1-l/2:2*dy:y1+l/2;
        [xx_s,yy_s] = meshgrid(x_source, y_source);
        E_src = (xx_s - x1 + 1i*(yy_s - y1));
        J = computeSource(E_src, k0, 2*dx, 2*dy);
        xidx=[round((x_source-mySub.subx(1))/dx)+1 round((x1-mySub.subx(1))/dx)+1];
        yidx=[round((y1-mySub.suby(1))/dy)+1 round((y_source-mySub.suby(1))/dy)+1];
        src_index = xidx + (yidx - 1) * mySub.Totalsubnx;
        mySub.b(src_index) = J(:);
    end
    if x2 > mySub.subx(1) && x2 < mySub.subx(2) && y2 > mySub.suby(1) && y2 < mySub.suby(2)
        x_source = x2-l/2:2*dx:x2+l/2;
        y_source = y2-l/2:2*dy:y2+l/2;
        [xx_s,yy_s] = meshgrid(x_source, y_source);
        E_src = (-(xx_s - x2) + 1i*(yy_s - y2));
        J = computeSource(E_src, k0, 2*dx, 2*dy);
        xidx=[round((x_source-mySub.subx(1))/dx)+1 round((x2-mySub.subx(1))/dx)+1];
        yidx=[round((y2-mySub.suby(1))/dy)+1 round((y_source-mySub.suby(1))/dy)+1];
        src_index = xidx + (yidx - 1) * mySub.Totalsubnx;
        mySub.b(src_index) = J(:);
    end
    fprintf('Worker %d: Complete source assembly\n', labindex);

    labBarrier;
    % Iterative solution
    Niter = 4;
    err = 0.01;
    for iter = 1:Niter
        if iter == 1
            Ez_new = mySub.Q * (mySub.U \ (mySub.L \ (mySub.P * mySub.b)));
        else
            Ez_new = mySub.Q * (mySub.U \ (mySub.L \ (mySub.P * mySub.Es)));
        end
        mySub.Ez = mySub.Ez + Ez_new;
        mySub.Es(:) = 0;
        % 1st synchronization: Ensure that all workers have completed this round of calculation
        fprintf('Worker %d: Complete the %d round of calculation, wait for other workers...\n', labindex, iter);
        labBarrier;

        % Step 2: Check convergence
        max_temp = max(abs(Ez_new));
        global_max = gop(@max, max_temp);  % Each lab gets the same global_max
        local_converged = global_max < err;  % Does the current lab consider it to have converged
        converged = gop(@or, local_converged);  % Any lab converges if it meets the conditions

        % If convergence occurs, all workers exit the iteration
        if converged
            break;
        end

        %Source transfer
        dir_list = {'left', 'right', 'top', 'bottom'};
        for d = 1:numel(dir_list)
            dir = dir_list{d};
            if isfield(mySub.neighbors, dir)
                nbr = mySub.neighbors.(dir);
                send_data = Ez_new(nbr.send);
                nbrID = nbr.id;
                switch dir
                    case {'left','right'}, parity = mod(myID,2);
                    case {'top','bottom'}, rowID = floor((myID-1)/NbrSubdomainx); parity = mod(rowID,2);
                end
                if parity == 0
                    recv_data = labReceive(nbrID);
                    labSend(send_data, nbrID);
                else
                    labSend(send_data, nbrID);
                    recv_data = labReceive(nbrID);
                end
                mySub.Es(nbr.recv) = mySub.Es(nbr.recv) - mySub.neighbors.(dir).T * recv_data;
            end
        end
        %  2nd synchronization: Ensure that all workers complete source transfer
        fprintf('Worker %d: Complete the %d round of source transfer, wait for other workers...\n', labindex, iter);
        labBarrier;
    end
    fprintf('Worker %d: Complete the iteration, %d steps\n', labindex, iter);

    % Collect global solutions
    local_result = struct();
    local_result.id = mySub.id;
    local_result.Ez = mySub.Ez;
    local_result.global_index = mySub.global_index;
end
tspmd=toc;
fprintf('Parallel time:%f',tspmd);
%% Postprocess
Ez_global = zeros(Totalnx * Totalny, 1);

% Collect results from all workers
for worker_idx = 1:length(local_result)
    worker_result = local_result{worker_idx};
    if ~isempty(worker_result.global_index)
        Ez_global(worker_result.global_index) = Ez_global(worker_result.global_index) + worker_result.Ez;
    end
end

%% Plot
Ez2d=reshape(Ez_global,Totalnx,Totalny);
n = 256; 
r = [(0:n/2)'/(n/2); ones(n/2+1, 1)];
g = [(0:n/2)'/(n/2); (n/2:-1:0)'/(n/2)];
b = [ones(n/2+1, 1); (n/2:-1:0)'/(n/2)];

cmap = [r, g, b];
figure;hold on
Plota=Totala/scale;
Plotb=Totalb/scale;
Plotxmin=((Py-1)*0.5*a0+R)/scale;
Plotymin=0;
Plotxmax=Plotxmin+Plota;
Plotymax=Plotb;

imagesc(linspace(Plotxmin,Plotxmax,Totalnx), linspace(Plotymin,Plotymax,Totalny), real(Ez2d.'));
axis xy;
axis equal tight;
set(gca, 'YDir', 'normal');  
colorbar;
colormap(cmap)
axis off;
model = mphload('Waveguide.mph');
hold on
hgeom = mphgeom(model, 'geom1', ...
    'facelabels','off', ...
    'edgecolor','k', ...
    'facealpha',0);  

set(hgeom,'LineWidth',0.5);
uistack(hgeom,'top');  
caxis([-2e-7 2e-7]);
hold on
xline((interfacex(2)+PMLnx*dx)/scale+Plotxmin);
%% functions
function J = computeSource(E_src, k0, dx, dy)
J_left = (-E_src(1,1) + 2*E_src(2,1) - E_src(3,1)) / dy^2 - k0^2 * E_src(2,1);
J_right = (-E_src(1,3) + 2*E_src(2,3) - E_src(3,3)) / dy^2 - k0^2 * E_src(2,3);
J_bottom = (-E_src(1,1) + 2*E_src(1,2) - E_src(1,3)) / dx^2 - k0^2 * E_src(1,2);
J_top = (-E_src(3,1) + 2*E_src(3,2) - E_src(3,3)) / dx^2 - k0^2 * E_src(3,2);
J = [J_left; J_bottom; J_right; J_top];
end
