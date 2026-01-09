%Set the parameters of the subdomain data structure
function [Subdomain] = SubParameters(Subdomain,k0,Totalnx,Totalny,epsilon,dx,dy)
PMLnx=8;
PMLny=8;
subx=Subdomain.subx;
suby=Subdomain.suby;
Totalsubnx= Subdomain.Totalsubnx;
Totalsubny= Subdomain.Totalsubny;
Subdomain.global_index=zeros(1,Totalsubnx*Totalsubny);
% Subdomain.Esmap=zeros(1,9*(Totalsubny));
% Subdomain.Efmap=zeros(1,9*Totalsubny);
Subdomain.eps=zeros(Totalsubnx,Totalsubny);
%% Get global_index and eps
for jj=1:Totalsubny
    for ii=1:Totalsubnx
        Subdomain.global_index(ii+(jj-1)*Totalsubnx)=ii+round(subx(1)/dx)+round(suby(1)/dy)*Totalnx+(jj-1)*Totalnx;
        Subdomain.eps(ii,jj)=epsilon(ii+round(subx(1)/dx),jj+round(suby(1)/dy));
    end
end
%% Get send and recv index of neighbors
% row first col sencond
for jj=1:Totalsubny
    for ii=1:9
        if isfield(Subdomain.neighbors,'left')
            Sleftindex=1:PMLnx+1;
            Rleftindex=PMLnx+1:PMLnx+1+PMLnx;
            Subdomain.neighbors.left.send(ii+(jj-1)*9)=Sleftindex(ii)+(jj-1)*Totalsubnx;
            Subdomain.neighbors.left.recv(ii+(jj-1)*9)=Rleftindex(ii)+(jj-1)*Totalsubnx;
        end
        if isfield(Subdomain.neighbors,'right')
            Srightindex=Totalsubnx-PMLnx:Totalsubnx;
            Rrightindex=Totalsubnx-2*PMLnx:Totalsubnx-PMLnx;
            Subdomain.neighbors.right.send(ii+(jj-1)*9)=Srightindex(ii)+(jj-1)*Totalsubnx;
            Subdomain.neighbors.right.recv(ii+(jj-1)*9)=Rrightindex(ii)+(jj-1)*Totalsubnx;
        end
    end
end

for jj=1:9
    for ii=1:Totalsubnx
        if isfield(Subdomain.neighbors,'top')
            Stopindex=Totalsubny-PMLny:Totalsubny;
            Rtopindex=Totalsubny-2*PMLny:Totalsubny-PMLny;
            Subdomain.neighbors.top.send(ii+(jj-1)*Totalsubnx)=ii+(Stopindex(jj)-1)*Totalsubnx;
            Subdomain.neighbors.top.recv(ii+(jj-1)*Totalsubnx)=ii+(Rtopindex(jj)-1)*Totalsubnx;
        end
        if isfield(Subdomain.neighbors,'bottom')
            Sbottomindex=1:PMLny+1;
            Rbottomindex=PMLny+1:PMLny+1+PMLny;
            Subdomain.neighbors.bottom.send(ii+(jj-1)*Totalsubnx)=ii+(Sbottomindex(jj)-1)*Totalsubnx;
            Subdomain.neighbors.bottom.recv(ii+(jj-1)*Totalsubnx)=ii+(Rbottomindex(jj)-1)*Totalsubnx;
        end
    end
end

%% T matrix
if isfield(Subdomain.neighbors,'left')||isfield(Subdomain.neighbors,'right')
    I=zeros(5,(9*Totalsubny));
    J=zeros(5,(9*Totalsubny));
    K=zeros(5,(9*Totalsubny));
    %  row first col second
    for n=PMLny+1:Totalsubny-PMLny;
        for m=2:9-1;
            index=(n-1)*9+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index+1,index-1,index+9,index-9];
            K(:,index)=[(-2/dx^2)+(-2/dy^2),1/dx^2,1/dx^2,...
                1/dy^2,1/dy^2];
        end
    end
    for m=[1,9]
        for n=1:Totalsubny
            index=(n-1)*9+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index,index,index,index];
            K(:,index)=[0,0,0,0,0];
        end
    end
    for n=[1:PMLny,Totalsubny-PMLny+1:Totalsubny]
        for m=1:9
            index=(n-1)*9+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index,index,index,index];
            K(:,index)=[0,0,0,0,0];
        end
    end
    if isfield(Subdomain.neighbors,'left')
        eps=zeros(9,Totalsubny);
        eps(2:end-1,PMLny+1:Totalsubny-PMLny)=Subdomain.eps(PMLnx+2:PMLnx*2,PMLny+1:Totalsubny-PMLny);
        epsvector=reshape(eps,9*Totalsubny,1);
        Subdomain.neighbors.left.T=sparse(I,J,K)+spdiags(epsvector*k0^2, 0, 9*Totalsubny,9*Totalsubny);
    end
    if isfield(Subdomain.neighbors,'right')
        eps=zeros(9,Totalsubny);
        eps(2:end-1,PMLny+1:Totalsubny-PMLny)=Subdomain.eps(Totalsubnx-PMLnx*2+1:Totalsubnx-PMLnx-1,PMLny+1:Totalsubny-PMLny);
        epsvector=reshape(eps,9*Totalsubny,1);
        Subdomain.neighbors.right.T=sparse(I,J,K)+spdiags(epsvector*k0^2, 0, 9*Totalsubny,9*Totalsubny);
    end
end

if isfield(Subdomain.neighbors,'top')||isfield(Subdomain.neighbors,'bottom')
    I=zeros(5,(9*Totalsubnx));
    J=zeros(5,(9*Totalsubnx));
    K=zeros(5,(9*Totalsubnx));
    %  row first col second
    for n=2:9-1
        for m=PMLnx+1:Totalsubnx-PMLnx
            index=(n-1)*Totalsubnx+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index+1,index-1,index+Totalsubnx,index-Totalsubnx];
            K(:,index)=[(-2/dx^2)+(-2/dy^2),1/dx^2,1/dx^2,...
                1/dy^2,1/dy^2];
        end
    end
    for n=[1,9]
        for m=1:Totalsubnx
            index=(n-1)*Totalsubnx+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index,index,index,index];
            K(:,index)=[0,0,0,0,0];
        end
    end
    for m=[1:PMLnx,Totalsubnx-PMLnx+1:Totalsubnx]
        for n=1:9
            index=(n-1)*Totalsubnx+m;
            I(:,index)=[index,index,index,index,index];
            J(:,index)=[index,index,index,index,index];
            K(:,index)=[0,0,0,0,0];
        end
    end
    if isfield(Subdomain.neighbors,'bottom')
        eps=zeros(Totalsubnx,9);
        eps(PMLnx+1:Totalsubnx-PMLnx,2:end-1)=Subdomain.eps(PMLnx+1:Totalsubnx-PMLnx,PMLny+2:PMLny*2);
        epsvector=reshape(eps,9*Totalsubnx,1);
        Subdomain.neighbors.bottom.T=sparse(I,J,K)+spdiags(epsvector*k0^2, 0, 9*Totalsubnx,9*Totalsubnx);
    end
    if isfield(Subdomain.neighbors,'top')
        eps=zeros(Totalsubnx,9);
        eps(PMLnx+1:Totalsubnx-PMLnx,2:end-1)=Subdomain.eps(PMLnx+1:Totalsubnx-PMLnx,Totalsubny-PMLny*2+1:Totalsubny-PMLny-1);
        epsvector=reshape(eps,9*Totalsubnx,1);
        Subdomain.neighbors.top.T=sparse(I,J,K)+spdiags(epsvector*k0^2, 0, 9*Totalsubnx,9*Totalsubnx);
    end
end

end