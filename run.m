addpath('_Grid_Data')
addpath('Grid_and_Cell_Metrics')
%% Properties/Initial Conditions
fluid.R = 287;        % Ideal Gas Constant [J/kg*K]
fluid.cp = 1005;      % Specific Heat [J/kg*k]
fluid.gamma = 1.4;    % Specific Heat Ratio
fluid.cv = fluid.cp/fluid.gamma;

free_stream.P_ref = 11664;  % Static Pressure [Pa]
free_stream.T_ref = 216.7;  % Temperature [K]
free_stream.M_ref = 3;      % Mach Number
free_stream.u_ref = 885;    % Axial Velocity [m/s]
free_stream.c_ref = 295;    % Speed of Sound [m/s]
free_stream.rho_ref = free_stream.P_ref/(fluid.R*free_stream.T_ref);

L = .015;   % Reference Length Scale [m]

%% Grid Setup
grid = setup_grid('g321x033uf.dat', 'plotFigs', false, 'Lref', L);
grid = cell_volume(grid);
grid = cell_face_area(grid);

%% Setup Initial Conditions
rho_et = free_stream.P_ref/(fluid.gamma-1) +...
    0.5*free_stream.rho_ref*free_stream.u_ref^2;

Q.q1 = grid.deltaV .* free_stream.rho_ref .* ones(grid.nx+1, grid.ny+1);    % rho
Q.q2 = grid.deltaV .* free_stream.rho_ref .*...
    free_stream.u_ref .* ones(grid.nx+1, grid.ny+1);                        % rho*u
Q.q3 = grid.deltaV .* zeros(grid.nx+1, grid.ny+1);                          % rho*v
Q.q4 = grid.deltaV .* rho_et .* ones(grid.nx+1, grid.ny+1);                 % rho*et

%%
for i = 1:100
    deltaT = time_step(Q, grid, fluid, 0.5);
    Q = apply_BC(Q, grid, fluid);
    [E, F] = roe_fds(Q, grid, fluid);
    [Q_new, E_diff, F_diff] = flux_diff(Q, E, F, grid, deltaT);

    Q = Q_new;

end

%%
u = Q.q2 ./ Q.q1;
v = Q.q3 ./ Q.q1;

contourf(grid.xc, grid.yc, u/885.13)