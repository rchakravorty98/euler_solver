function grid = cell_face_interp(grid)
%%% Interpolate Cell Face Metrics to Cell Center (used in Time Step)
nx = grid.nx;
ny = grid.ny;

%% xi direction
f = scatteredInterpolant(grid.xi.x(:), grid.xi.y(:), grid.xi.Sx(:));
grid.xi.Sx_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));

f = scatteredInterpolant(grid.xi.x(:), grid.xi.y(:), grid.xi.Sy(:));
grid.xi.Sy_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));

f = scatteredInterpolant(grid.xi.x(:), grid.xi.y(:), grid.xi.S(:));
grid.xi.S_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));

%% eta direction
f = scatteredInterpolant(grid.eta.x(:), grid.eta.y(:), grid.eta.Sx(:));
grid.eta.Sx_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));

f = scatteredInterpolant(grid.eta.x(:), grid.eta.y(:), grid.eta.Sy(:));
grid.eta.Sy_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));

f = scatteredInterpolant(grid.eta.x(:), grid.eta.y(:), grid.eta.S(:));
grid.eta.S_avg = f(grid.xc(2:nx,2:ny), grid.yc(2:nx,2:ny));


end