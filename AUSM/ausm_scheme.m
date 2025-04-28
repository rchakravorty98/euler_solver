function [E, F] = ausm_scheme(Q, grid, fluid, eps, k)
nx = grid.nx;
ny = grid.ny;

%% xi direction
i_idx = 2:nx+1;
j_idx = 2:ny;   % Need fluxes over primary grid

[Q_l, Q_r] = muscl(Q, grid, i_idx, j_idx, eps, k, 'xi');

n_x = grid.xi.Sx./grid.xi.S;
n_y = grid.xi.Sy./grid.xi.S;
deltaV_l = grid.deltaV(i_idx-1,j_idx);
deltaV_r = grid.deltaV(i_idx,j_idx);

[E.e1, E.e2, E.e3, E.e4] = ausm_flux_construction(Q_l, Q_r, deltaV_l, deltaV_r, n_x, n_y, fluid, 'xi');

%% eta direction
i_idx = 2:nx;   % Need fluxes over primary grid
j_idx = 2:ny+1;

[Q_l, Q_r] = muscl(Q, grid, i_idx, j_idx, eps, k, 'eta');

n_x = grid.eta.Sx./grid.eta.S;
n_y = grid.eta.Sy./grid.eta.S;
deltaV_l = grid.deltaV(i_idx,j_idx-1);
deltaV_r = grid.deltaV(i_idx,j_idx);

[F.f1, F.f2, F.f3, F.f4] = ausm_flux_construction(Q_l, Q_r, deltaV_l, deltaV_r, n_x, n_y, fluid, 'eta');

end