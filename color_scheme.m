function cmap = color_scheme(n)

% basic colors from yellow to dark red
cmap_ori = [1, 0.85, 0.46;
            1, 0.70, 0.30;
            1, 0.55, 0.24;
            0.99, 0.30, 0.17;
            0.89, 0.10,	0.11;
            0.70, 0, 0.15];

% linearly interpolate
cmap1 = interp1(1:6,cmap_ori(:,1)',linspace(1,6,n))';
cmap2 = interp1(1:6,cmap_ori(:,2)',linspace(1,6,n))';
cmap3 = interp1(1:6,cmap_ori(:,3)',linspace(1,6,n))';

cmap = [cmap1, cmap2, cmap3];
