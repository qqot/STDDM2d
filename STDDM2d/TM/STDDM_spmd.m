%  A source transfer domain decomposition method to calculate TM scattering
%  problem of 2d dielectric cylinder
%  The code was written by Zhanwen Wang from Zhejiang University
%  Email: zw_wang@zju.edu.cn

clear;
clc;

%% Parameters
lambda0=1;
PMLnx=8;               % Number of PML layers
PMLny=8;               % Number of PML layers
dx=0.01*lambda0;       % Length of the mesh
dy=0.01*lambda0;       % Width of the mesh
a = 80*lambda0;        %Length
b = 80*lambda0;
nx=round(a/dx+1);      % Number of grids (x-direction) in the region of interest
ny=round(b/dy+1);      % Number of grids (y-direction) in the region of interest

c=3e8;                 % Speed of light


k0=2*pi/lambda0;       % Wavenumber
omega0=k0*c;           % Angular frequency
eps0=8.854e-12;
eps_1=1;
eps_2=4;

Totalnx=nx+2*PMLnx;
Totalny=ny+2*PMLny;
Totala=a+2*PMLnx*dx;
Totalb=b+2*PMLny*dy;

r=9*lambda0;  %radius of dielectric cylinder
x0=a/2+PMLnx*dx;
y0=b/2+PMLny*dy;
%epsilon
epsilon_center=Geteps(Totalnx,Totalny,x0,y0,r,dx,dy,eps_1, eps_2);
epsilon=ones(Totalnx,Totalny);
i_range = 2:Totalnx-1;
j_range = 2:Totalny-1;

% Average the four surrounding cell centers for each grid point
epsilon(i_range, j_range) = (epsilon_center(i_range-1, j_range-1) + ...
    epsilon_center(i_range-1, j_range) + ...
    epsilon_center(i_range, j_range-1) + ...
    epsilon_center(i_range, j_range)) / 4;

%% Domain Decomposition
interfacex=roundn(linspace(0,a,5),-1);
interfacey=roundn(linspace(0,b,5),-1);
NbrSubdomainx=length(interfacex)-1;NbrSubdomainy=length(interfacey)-1;
Ndomain=NbrSubdomainx*NbrSubdomainy;
Subdomains=cell(Ndomain,1);
idx_matrix = reshape(1:NbrSubdomainx*NbrSubdomainy, NbrSubdomainx, NbrSubdomainy);

% ---------- Start the parallel pool ----------
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
    mySub = AssemblyofEqu(mySub, k0, dx, dy,dz);
    mySub.Es=zeros(mySub.Totalsubnx*mySub.Totalsubny,1);
    mySub.Ez=zeros(mySub.Totalsubnx*mySub.Totalsubny,1);
    fprintf('Worker %d: Complete the matrix assembly\n', labindex);
    % Sparse LU factorization
    [mySub.L, mySub.U, mySub.P, mySub.Q] = lu(mySub.A);
    fprintf('Worker %d: Complete factorization\n', labindex);

    labBarrier;
    % 迭代求解
    Niter = 30;
    err = 1;
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

        %         % 第二步：检查收敛性（仅在非第一次迭代时）
        %         max_temp = max(abs(Ez_new));
        %         global_max = gop(@max, max_temp);  % 每个lab都获得相同global_max
        %         local_converged = global_max < err;  % 当前lab是否认为已收敛
        %         converged = gop(@or, local_converged);  % 任一lab满足条件即收敛
        %
        %         % 如果收敛，所有worker都退出迭代
        %         if converged
        %             break;
        %         end


        % Source transfer
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
        % 2nd synchronization: Ensure that all workers complete source transfer
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
fprintf('Parallel time: %f',tspmd);
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
figure;
imagesc(linspace(0,Totala,Totalnx), linspace(0,Totalb,Totalny), abs(Ez2d.'));
axis equal;
axis xy;
colorbar;
colormap jet
%% RCS
[RCS] = RCS_2d(Ez2d,lambda0,Totalnx,Totalny,dx,dy);

figure
plot(linspace(0, 360, 361),RCS,'--','LineWidth',1.5);
grid on;
