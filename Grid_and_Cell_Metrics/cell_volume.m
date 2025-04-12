function grid = cell_volume(grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds cell volumes of the inputs grid
%
% Inputs:
%   - grid: Structure with following fields:
%
%       - x:    x grid points corresponding to cell vertices
%       - y:    y grid points corresponding to cell vertices
%       - nx:   Number of primary x grid vertices
%       - ny:   Number of primary y grid vertices
%
% Outputs (added as additional fields to grid input structure):
%   - xc:       x grid points corresponding to cell centroids
%   - yc:       y grid points corresponding to cell centroids
%   - deltaV:   Cell Volumes
%
% Since MATLAB indexing starts at 1:
%
%   - Grid Array: 1 ≤ i ≤ nx+2 & 1 ≤ j ≤ ny+2
%   - Primary Grid: 2 ≤ i ≤ nx+1 & 2 ≤ j ≤ ny+1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = grid.x;
y = grid.y;
nx = grid.nx;
ny = grid.ny;

deltaV = zeros(nx+1, ny+1);
xc = zeros(nx+1, ny+1); 
yc = zeros(nx+1, ny+1);
for i = 1:nx+1
    for j = 1:ny+1
        % i-1/2, j-1/2
        x1 = x(i,j);
        y1 = y(i,j);

        % i+1/2, j-1/2
        x2 = x(i+1, j);
        y2 = y(i+1, j);

        % i+1/2, j+1/2
        x3 = x(i+1, j+1);
        y3 = y(i+1, j+1);

        % i-1/2, j+1/2
        x4 = x(i,j+1);
        y4 = y(i,j+1);

        r1 = [(x3-x1) (y3-y1)];
        r2 = [(x4-x2) (y4-y2)];

        deltaV(i,j) = 0.5*abs(r1(1)*r2(2) - r2(1)*r1(2));

        %%% Find Centroid of Cell (i.e. Centroid of Arbitrary
        %%% Quadrilateral)
        xc1 = mean([x1 x2 x3]);
        yc1 = mean([y1 y2 y3]);
        A1 = 0.5*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));

        xc2 = mean([x1 x3 x4]);
        yc2 = mean([y1 y3 y4]);
        A2 = 0.5*(x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3));

        xc(i,j) = (xc1*A1+xc2*A2)/(A1+A2);
        yc(i,j) = (yc1*A1+yc2*A2)/(A1+A2);

    end
end

grid.xc = xc;
grid.yc = yc;
grid.deltaV = deltaV;

end