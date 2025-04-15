function deltaT = time_step(Q, grid, fluid, CFL_max)
nx = grid.nx;
ny = grid.ny;

u = Q.q2 ./ Q.q1;
v = Q.q3 ./ Q.q1;
et = Q.q4 ./ Q.q1;
T = (et - 0.5*(u.^2+v.^2))/fluid.cv;
c = (fluid.gamma * fluid.R * T).^(1/2);

%% xi Direction
n_x = grid.xi.Sx_avg ./ grid.xi.S_avg;
n_y = grid.xi.Sy_avg ./ grid.xi.S_avg;

contra_vel = u(2:nx, 2:ny).*n_x + v(2:nx, 2:ny).*n_y;     % Contravariant Velocity

xi_x = grid.xi.Sx_avg ./ grid.deltaV(2:nx, 2:ny);
xi_y = grid.xi.Sy_avg ./ grid.deltaV(2:nx, 2:ny);
C1 = CFL_max*(1./((abs(contra_vel) + c(2:nx, 2:ny)) .* (xi_x.^2 + xi_y.^2).^(1/2)));

%% eta Direction
n_x = grid.eta.Sx_avg ./ grid.eta.S_avg;
n_y = grid.eta.Sy_avg ./ grid.eta.S_avg;

contra_vel = u(2:nx, 2:ny).*n_x + v(2:nx, 2:ny).*n_y;     % Contravariant Velocity

eta_x = grid.eta.Sx_avg ./ grid.deltaV(2:nx, 2:ny);
eta_y = grid.eta.Sy_avg ./ grid.deltaV(2:nx, 2:ny);
C2 = CFL_max*(1./((abs(contra_vel) + c(2:nx, 2:ny)) .* (eta_x.^2 + eta_y.^2).^(1/2)));

%% Time Step
% deltaT = min([min(C1) min(C2)]);
deltaT = min(C1, C2);

end