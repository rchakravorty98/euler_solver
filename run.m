addpath('Grid_Data')
addpath('Grid_and_Cell_Metrics')
addpath('AUSM')
addpath('General')

clear;close all;
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
grid = setup_grid("g641x065uf.dat", 'plotFigs', false, 'Lref', L);
grid = cell_volume(grid);
grid = cell_face_area(grid);
grid = cell_face_interp(grid);

%% Setup Initial Conditions
free_stream.rho_et = free_stream.P_ref/(fluid.gamma-1) +...
    0.5*free_stream.rho_ref*free_stream.u_ref^2;

Q.q1 = grid.deltaV .* free_stream.rho_ref .* ones(grid.nx+1, grid.ny+1);    % rho
Q.q2 = grid.deltaV .* free_stream.rho_ref .*...
    free_stream.u_ref .* ones(grid.nx+1, grid.ny+1);                        % rho*u
Q.q3 = grid.deltaV .* zeros(grid.nx+1, grid.ny+1);                          % rho*v
Q.q4 = grid.deltaV .* free_stream.rho_et .* ones(grid.nx+1, grid.ny+1);     % rho*et

%% Main Loop
max_iter = 100000;
res_history = zeros(1, max_iter);
eps = 1;
for i = 1:max_iter
    
    % if i >= 100
        CFLmax = 0.05;
        eps = 1;
        k = 0.5;
    % else
        % CFLmax = 0.3;
        % eps = 0;
        % k = 0;
    % end


    Q = apply_BC(Q, grid, fluid);
    deltaT = time_step(Q, grid, fluid, CFLmax);

    [E, F] = ausm_scheme(Q, grid, fluid, eps, k);
    [Q_new, E_diff, F_diff] = flux_diff(Q, E, F, grid, deltaT);
    
    [res_history(i), grid_error] = residual(Q, Q_new, grid, free_stream);

    Q = Q_new;

    fprintf('Iteration Num: %d | Worst: %d \n ', i, res_history(i))
    
    if any(~isreal(deltaT))
        break
    end

end

%%
[rho, u, v, et, P, T] = Q_to_primitive(Q.q1, Q.q2, Q.q3, Q.q4, grid.deltaV, fluid);
c = (fluid.gamma .* fluid.R .* T).^(1/2);

figure()
semilogy(1:numel(res_history), res_history, 'LineWidth', 2)
grid on
xlabel('Iteration Number')
ylabel('L_{\infty}-Norm')

fig = figure(2);
tiledlayout(4,1)

nexttile
u_norm = u/free_stream.u_ref;
contourf(grid.xc(2:grid.nx, 2:grid.ny), grid.yc(2:grid.nx, 2:grid.ny), u_norm(2:grid.nx, 2:grid.ny),...
    [min(min(u_norm)):0.0005:max(max(u_norm))], 'LineColor', 'none')
colormap turbo
colorbar()
title('U/U_{ref}')

nexttile
p_norm = (P-free_stream.P_ref) ./ (free_stream.rho_ref*free_stream.u_ref^2);
contourf(grid.xc(2:grid.nx, 2:grid.ny), grid.yc(2:grid.nx, 2:grid.ny), p_norm(2:grid.nx, 2:grid.ny),...
    [min(min(p_norm)):0.0005:max(max(p_norm))], 'LineColor', 'none')
colormap turbo
colorbar()
title('(P-P_{ref})/(\rho_{ref}*U_{ref}^2)')

nexttile
T_norm = (T-free_stream.T_ref) ./ ((free_stream.c_ref.^2)/1004.5);
contourf(grid.xc(2:grid.nx, 2:grid.ny), grid.yc(2:grid.nx, 2:grid.ny), T_norm(2:grid.nx, 2:grid.ny),...
    [min(min(T_norm)):0.0005:max(max(T_norm))], 'LineColor', 'none')
colormap turbo
colorbar()
title('(T-T_{ref})/(c_{ref}^{2}*cp_{ref})')

nexttile
rho_norm = rho ./ free_stream.rho_ref;
contourf(grid.xc(2:grid.nx, 2:grid.ny), grid.yc(2:grid.nx, 2:grid.ny), rho_norm(2:grid.nx, 2:grid.ny),...
    [min(min(rho_norm)):0.0005:max(max(rho_norm))], 'LineColor', 'none')
colormap turbo
colorbar()
title('\rho/\rho_{ref}')
fig.Position = [0 0 fig.Position(3)*3.25 fig.Position(4)*3];
