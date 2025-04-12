function Q = apply_BC(Q, grid, fluid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applies Boundary Conditions to Q Matrix. BC's applied:
%
%   - Slip Wall
%   - Adiabatic Wall
%   - 0 Pressure Gradient at Wall
%   - Q at nx+1 (halo cell) = nx 
%
% Inputs:
%   - Q:        Q vector at every cell
%   - grid:     Structure containing all grid parameters
%   - fluid:    Structure containing fluid properties
%
% Outputs (added as additional fields to grid input structure):
%   - Q:    Updated Q with BC's applied
%
% Since MATLAB indexing starts at 1:
%
%   - Grid Array: 1 ≤ i ≤ nx+2 & 1 ≤ j ≤ ny+2
%   - Primary Grid: 2 ≤ i ≤ nx+1 & 2 ≤ j ≤ ny+1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = grid.nx;
ny = grid.ny;

cv = fluid.cv;
R = fluid.R;
gamma = fluid.gamma;

%% Slip Wall
S_eta_x = grid.eta.Sx;
S_eta_y = grid.eta.Sy;

% Wall Indices:
%
%   - Top Wall (2 ≤ i ≤ nx at j = ny+1)
%   - Bottom Wall (2 ≤ i ≤ nx at j = 2) 
%
% Note: Cell face areas only span primary grid
for j = [2 ny+1] 

    if j == 2
        int_ind = j;
        halo_ind = j-1;
    else
        int_ind = j-1;
        halo_ind = j;
    end
    
    C1 = 1 ./ (S_eta_x(:,j-1).^2 + S_eta_y(:,j-1).^2);
    C2 = S_eta_x(:,j-1).^2 - S_eta_y(:,j-1).^2;
    C3 = S_eta_x(:,j-1) .* S_eta_y(:,j-1);

    % Interior Velocities
    u_int = Q.q2(2:nx, int_ind) ./ Q.q1(2:nx, int_ind);
    v_int = Q.q3(2:nx, int_ind) ./ Q.q1(2:nx, int_ind);

    % Required Velocity Components in Halo Cells
    u_o = C1 .* (-u_int .* C2 -  2 .* v_int .* C3);
    v_o = C1 .* (v_int .* C2 -  2 .* u_int .* C3);

    % Temperature and Pressure in Interior Cell
    rho_int = Q.q1(2:nx, int_ind) ./ grid.deltaV(2:nx, int_ind);    % Interior Density
    et = Q.q4(2:nx, int_ind) ./ Q.q1(2:nx, int_ind);                % Interior Total Energy
    T_int = (et - 0.5*(u_int.^2+v_int.^2))/cv;                      % Interior Temperature
    P_int = rho_int .* R .* T_int;                                  % Interior Pressure

    % Halo Cell
    rho_o = P_int ./ (R .* T_int);
    rho_et_o = P_int ./ (gamma-1) + 0.5 .* rho_o .* (u_o.^2 + v_o.^2);

    % Replace Halo Cell
    Q.q1(2:nx, halo_ind) = rho_o .* grid.deltaV(2:nx, halo_ind);
    Q.q2(2:nx, halo_ind) = rho_o .* u_o .* grid.deltaV(2:nx, halo_ind);
    Q.q3(2:nx, halo_ind) = rho_o .* v_o .* grid.deltaV(2:nx, halo_ind);
    Q.q4(2:nx, halo_ind) = rho_et_o .* grid.deltaV(2:nx, halo_ind);
end

%% Exit Plane
Q.q1(nx+1, :) = grid.deltaV(nx+1,:) .* (Q.q1(nx,:) ./ grid.deltaV(nx,:));
Q.q2(nx+1, :) = grid.deltaV(nx+1,:) .* (Q.q2(nx,:) ./ grid.deltaV(nx,:));
Q.q3(nx+1, :) = grid.deltaV(nx+1,:) .* (Q.q3(nx,:) ./ grid.deltaV(nx,:));
Q.q4(nx+1, :) = grid.deltaV(nx+1,:) .* (Q.q4(nx,:) ./ grid.deltaV(nx,:));

end