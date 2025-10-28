% subcell integration method to accurately compute the effective relative permittivity
function eps_r = Geteps(Totalnx,Totalny,cx,cy,radius,dx,dy,eps_outside, eps_inside)
% Grid center coordinates
x = ((1:Totalnx-1) - 0.5) * dx;
y = ((1:Totalny-1) - 0.5) * dy;
% The coordinate matrix of the main grid center
[X, Y] = meshgrid(x, y);
X=X.';Y=Y.';
    % Initialize to the background relative electric permittivity
    eps_r = eps_outside * ones(Totalnx-1, Totalny-1);
    % Calculate the coordinate offset of the four vertices of each cell in the main grid
    half_dx = dx / 2;
    half_dy = dy / 2;
    % Four vertex coordinates
    x1 = X - half_dx; y1 = Y - half_dy;
    x2 = X + half_dx; y2 = Y - half_dy;
    x3 = X - half_dx; y3 = Y + half_dy;
    x4 = X + half_dx; y4 = Y + half_dy;
    % The distances from the four vertices to the center of the circle
    r1 = sqrt((x1 - cx).^2 + (y1 - cy).^2);
    r2 = sqrt((x2 - cx).^2 + (y2 - cy).^2);
    r3 = sqrt((x3 - cx).^2 + (y3 - cy).^2);
    r4 = sqrt((x4 - cx).^2 + (y4 - cy).^2);
    % Determine whether the unit is completely inside or outside the circle
    all_outside = (r1 > radius) & (r2 > radius) & (r3 > radius) & (r4 > radius);
    all_inside  = (r1 < radius) & (r2 < radius) & (r3 < radius) & (r4 < radius);
    % Set all the units directly within the circle
    eps_r(all_inside) = eps_inside;

    need_subgrid = ~(all_inside | all_outside);
    % Extract the index that needs to be processed
    [ix_list, iy_list] = find(need_subgrid);

    for nn=1:length(ix_list)
        ix=ix_list(nn);
        iy=iy_list(nn);
        for ii=-10:10
            for jj=-10:10
                xdist=(ix-0.5-ii*0.2)*dx-cx;
                ydist=(iy-0.5-jj*0.2)*dy-cy;
                dist=sqrt(xdist.^2+ydist.^2);
                if dist<=radius
                    eps_r(ix,iy)=eps_r(ix,iy)+1/441*(eps_inside-eps_outside);
                end
            end
        end
    end

end