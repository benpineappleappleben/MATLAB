function breakaway_cube(S, n, L, figure_number)
    clim = [0, max(S(:))];
    cmap = crameri('batlow');
    figure(figure_number); clf; hold on;
    colormap(cmap);
    axis([0,L(1),0,L(2),0,L(3)]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Breakaway Cube Visualization');
    [nx, ny, nz] = size(S);

    x = linspace(0, L(1), nx);
    y = linspace(0, L(2), ny);

    v = zeros(n,1);
    for i47 = 0:n-1
        if i47 == 0
            v(i47+1) = 1;
        else
            v(i47+1) = 1 + round(i47*(nz-1)/(n-1));
        end
    end

    for k1 = 1:n
        slice_data = S(:,:,v(k1));
        zpos_scaled = (v(k1)-1)/(nz-1) * L(3);
        surf(x, y, zpos_scaled*ones(nx,ny), slice_data','EdgeColor','none','FaceColor','interp');
    end

    colorbar;
    caxis(clim);
    view(15,5);
    camlight left; lighting gouraud;
    drawnow;
end