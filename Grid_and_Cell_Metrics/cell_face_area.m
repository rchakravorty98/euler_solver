function grid = cell_face_area(grid)
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
%   - xi/eta:   Contain the following fields for the xi and eta direction
%
%       - S:    Cell Face Area
%       - Sx:   Projected Cell Face Area in x
%       - Sy:   Projected Cell Face Area in y
%       - x:    x Location of Midpoint of Cell Faces Perpendicular to xi/eta
%       - y:    y Location of Midpoint of Cell Faces Perpendicular to xi/eta
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

%% xi Direction
S_xi = zeros(nx, ny-1);
S_xi_x = zeros(nx, ny-1);
S_xi_y = zeros(nx, ny-1);
x_xi = zeros(nx, ny-1);
y_xi = zeros(nx, ny-1);
for i = 2:nx+1
    for j = 2:ny 
        S_xi_x(i-1,j-1) = y(i,j+1)-y(i,j);
        S_xi_y(i-1,j-1) = -(x(i,j+1)-x(i,j));
        S_xi(i-1,j-1) = ((S_xi_x(i-1,j-1))^2 + (S_xi_y(i-1,j-1))^2)^(1/2);

        %%% Location of Midpoint of Cell Faces Perpendicular to xi Direction
        x_xi(i-1,j-1) = (x(i,j+1)+x(i,j))/2;
        y_xi(i-1,j-1) = (y(i,j+1)+y(i,j))/2;
    end
end

xi.S = S_xi;
xi.Sx = S_xi_x;
xi.Sy = S_xi_y;
xi.x = x_xi;
xi.y = y_xi;

%% eta Direction
S_eta = zeros(nx-1, ny);
S_eta_x = zeros(nx-1, ny);
S_eta_y = zeros(nx-1, ny);
x_eta = zeros(nx-1, ny);
y_eta = zeros(nx-1, ny);
for i = 2:nx
    for j = 2:ny+1 
        S_eta_x(i-1,j-1) = -(y(i+1,j)-y(i,j));
        S_eta_y(i-1,j-1) = x(i+1,j)-x(i,j);
        S_eta(i-1,j-1) = ((S_eta_x(i-1,j-1))^2 + (S_eta_y(i-1,j-1))^2)^(1/2);

        %%% Location of Midpoint of Cell Faces Perpendicular to eta Direction
        x_eta(i-1,j-1) = (x(i+1,j)+x(i,j))/2;
        y_eta(i-1,j-1) = (y(i+1,j)+y(i,j))/2;
    end
end

eta.S = S_eta;
eta.Sx = S_eta_x;
eta.Sy = S_eta_y;
eta.x = x_eta;
eta.y = y_eta;

%% Ouput
grid.xi = xi;
grid.eta = eta;

end