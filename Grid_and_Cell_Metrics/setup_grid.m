function grid = setup_grid(filename, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in ASCII file and outputs mesh with Halo Cells
%
% Inputs:
%   - filename: Name of file containing x/y pairs
%
% Outputs:
%   - grid: Structure with following fields
%
%       - x:    x grid points corresponding to cell vertices with Halo Cells
%       - y:    y grid points corresponding to cell vertices with Halo Cells
%       - nx:   Number of primary x grid vertices
%       - ny:   Number of primary y grid vertices
%
% Since MATLAB indexing starts at 1:
%
%   - Grid Array: 1 ≤ i ≤ nx+2 & 1 ≤ j ≤ ny+2
%   - Primary Grid: 2 ≤ i ≤ nx+1 & 2 ≤ j ≤ ny+1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    filename {mustBeMember(filename, {'test_grid.dat', 'g641x065uf.dat', 'g321x033uf.dat'})}
    options.plotFigs logical = false; 
    options.Lref = 1;   % Used to dimensionalize grid
end
Lref = options.Lref;

%% Read File
fid = fopen(filename, 'r');

%%% Extract nx and ny
l = fgetl(fid); % Skip First Line
nx_ny = regexp(l, 'i=(\d+),\s*j=(\d+)', 'tokens');

nx = str2double(nx_ny{1}{1});
ny = str2double(nx_ny{1}{2});

%%% Extract x,y Pairs
data = fscanf(fid, '%f,%f', [2 Inf]); % Read x,y pairs
fclose(fid);

%% Rearrange x,y Pairs into Mesh
x_ij = zeros(nx, ny);
y_ij = zeros(nx, ny);

counter = 1;
for j = 1:ny
    for i = 1:nx
        x_ij(i,j) = data(1,counter);
        y_ij(i,j) = data(2,counter);
        
        counter = counter + 1;
    end
end
x_ij = x_ij*Lref;
y_ij = y_ij*Lref;

%% Construct Halo Cells
x_ij_halo = zeros(nx+2, ny+2);
y_ij_halo = zeros(nx+2, ny+2);

%%% Fill in Primary Grid
x_ij_halo(2:nx+1, 2:ny+1) = x_ij;
y_ij_halo(2:nx+1, 2:ny+1) = y_ij;

%%% Left
x_ij_halo(1, 2:ny+1) = 2*x_ij_halo(2, 2:ny+1) - x_ij_halo(3, 2:ny+1);
y_ij_halo(1, 2:ny+1) = y_ij_halo(2, 2:ny+1);

%%% Right
x_ij_halo(nx+2, 2:ny+1) = 2*x_ij_halo(nx+1, 2:ny+1) - x_ij_halo(nx, 2:ny+1);
y_ij_halo(nx+2, 2:ny+1) = y_ij_halo(nx+1, 2:ny+1);

%%% Bottom
x_ij_halo(1:nx+2, 1) = x_ij_halo(1:nx+2, 2);
y_ij_halo(1:nx+2, 1) = 2*y_ij_halo(1:nx+2, 2) - y_ij_halo(1:nx+2, 3);

%%% Top
x_ij_halo(1:nx+2, ny+2) = x_ij_halo(1:nx+2, ny+1);
y_ij_halo(1:nx+2, ny+2) = 2*y_ij_halo(1:nx+2, ny+1) - y_ij_halo(1:nx+2, ny);

%% Output
grid.x = x_ij_halo;
grid.y = y_ij_halo;
grid.nx = nx;
grid.ny = ny;

%% Plot
if options.plotFigs
    fig = figure();
    hold on

    %%% Ghost Cells
    plot(x_ij_halo(:,1:2), y_ij_halo(:,1:2), 'b', x_ij_halo(:,1:2)', y_ij_halo(:,1:2)', 'b')
    plot(x_ij_halo(:,ny+1:ny+2), y_ij_halo(:,ny+1:ny+2), 'b', x_ij_halo(:,ny+1:ny+2)', y_ij_halo(:,ny+1:ny+2)', 'b')
    plot(x_ij_halo(1:2,:), y_ij_halo(1:2,:), 'b', x_ij_halo(1:2,:)', y_ij_halo(1:2,:)', 'b')
    plot(x_ij_halo(nx+1:nx+2,:), y_ij_halo(nx+1:nx+2,:), 'b', x_ij_halo(nx+1:nx+2,:)', y_ij_halo(nx+1:nx+2,:)', 'b')
    xlim([min(min(x_ij_halo)) max(max(x_ij_halo))])
    ylim([min(min(y_ij_halo)) max(max(y_ij_halo))])

    %%% Primary Grid
    plot(x_ij, y_ij, 'k', x_ij', y_ij', 'k')
    xlabel('x/L')
    ylabel('y/L')
    title('Augmented Grid')
    fig.Position = [0 0 fig.Position(3)*3.25 fig.Position(4)];
end

end