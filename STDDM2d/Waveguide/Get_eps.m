function [eps] = Get_eps(a,b,dx,dy,delta,scale)
%Get epsilon from comsol
model=mphload('Waveguide.mph');
nx = round(a/dx + 1);
ny =round(b/dy+ 1);
[ix_grid, iy_grid] = ndgrid(1:nx-1, 1:ny-1);
x_grid = (ix_grid - 1) * dx+dx/2+delta;
y_grid = (iy_grid - 1) * dy+dy/2;
coords = [x_grid(:)'; y_grid(:)']; % 转成 2×N 的形式
coords=coords/scale;
n_iso = mphinterp(model, 'ewfd.n_iso', 'coord', coords);
eps=reshape(n_iso.^2,nx-1,ny-1);
end