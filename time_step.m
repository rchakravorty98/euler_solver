function deltaT = time_step(Q, grid, fluid, CFL_max)
nx = grid.nx;
ny = grid.ny;

u = Q.q2 ./ Q.q1;
v = Q.q3 ./ Q.q1;
et = Q.q4 ./ Q.q1;
T = (et - 0.5*(u.^2+v.^2))/fluid.cv;
c = (fluid.gamma * fluid.R * T).^(1/2);

%% xi Direction
%%% Average Cell Face Metrics to Cell Center
Sx_avg = (grid.xi.Sx(1:end-1,:) + grid.xi.Sx(2:end,:))/2;
Sy_avg = (grid.xi.Sy(1:end-1,:) + grid.xi.Sy(2:end,:))/2;
S_avg = (grid.xi.S(1:end-1,:) + grid.xi.S(2:end,:))/2;

n_x = Sx_avg ./ S_avg;
n_y = Sy_avg ./ S_avg;

contra_vel = u(2:nx, 2:ny).*n_x + v(2:nx, 2:ny).*n_y;     % Contravariant Velocity

xi_x = Sx_avg ./ grid.deltaV(2:nx, 2:ny);
xi_y = Sy_avg ./ grid.deltaV(2:nx, 2:ny);
C1 = CFL_max*(1./((abs(contra_vel) + c(2:nx, 2:ny)) .* (xi_x.^2 + xi_y.^2)));

%% eta Direction
%%% Average Cell Face Metrics to Cell Center
Sx_avg = (grid.eta.Sx(:,1:end-1) + grid.eta.Sx(:,2:end))/2;
Sy_avg = (grid.eta.Sy(:,1:end-1) + grid.eta.Sy(:,2:end))/2;
S_avg = (grid.eta.S(:,1:end-1) + grid.eta.S(:,2:end))/2;

n_x = Sx_avg ./ S_avg;
n_y = Sy_avg ./ S_avg;

contra_vel = u(2:nx, 2:ny).*n_x + v(2:nx, 2:ny).*n_y;     % Contravariant Velocity

eta_x = Sx_avg ./ grid.deltaV(2:nx, 2:ny);
eta_y = Sy_avg ./ grid.deltaV(2:nx, 2:ny);
C2 = CFL_max*(1./((abs(contra_vel) + c(2:nx, 2:ny)) .* (eta_x.^2 + eta_y.^2)));

%% Time Step
deltaT = min([min(C1) min(C2)]);

end