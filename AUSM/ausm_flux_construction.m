function [f1, f2, f3, f4] = ausm_flux_construction(Q_l, Q_r, deltaV_l, deltaV_r, n_x, n_y, fluid, dir)
gamma = fluid.gamma;
R = fluid.R;

%%% Constants
beta = 1/8;
Ku = 3/4;
Kp = 1/4;
sigma = 1;

%% Left State
[rho_l, u_l, v_l, et_l, P_l, T_l, ht_l] = Q_to_primitive(Q_l.q1, Q_l.q2,...
    Q_l.q3, Q_l.q4, deltaV_l, fluid);

contra_vel_left = u_l.* n_x + v_l.*n_y; 

psi_l = [ones(1,1,numel(Q_l.q1));...
    reshape(u_l, 1, 1, []);...
    reshape(v_l, 1, 1, []);...
    reshape(ht_l, 1, 1, [])];

%% Right State 
[rho_r, u_r, v_r, et_r, P_r, T_r, ht_r] = Q_to_primitive(Q_r.q1, Q_r.q2,...
    Q_r.q3, Q_r.q4, deltaV_r, fluid);

contra_vel_right = u_r.* n_x + v_r.*n_y;

psi_r = [ones(1,1,numel(Q_r.q1));...
    reshape(u_r, 1, 1, []);...
    reshape(v_r, 1, 1, []);...
    reshape(ht_r, 1, 1, [])];

%% Half States
% c_half = 0.5 * (c_l + c_r);
rho_half =  0.5 * (rho_l + rho_r);

crit_cl = ( 2 * ht_l * (gamma-1)/(gamma+1) ).^(1/2);
crit_cr = ( 2 * ht_r * (gamma-1)/(gamma+1) ).^(1/2);
c_hat_l = (crit_cl.^2) ./ max(crit_cl, contra_vel_left);
c_hat_r = (crit_cr.^2) ./ max(crit_cr, -contra_vel_right);
c_half = min(c_hat_r, c_hat_l);

M_l = contra_vel_left ./ c_half;
M_r = contra_vel_right ./ c_half;
M_bar2 = 0.5*(M_l.^2 + M_r.^2);

Mo = (min(1,max(M_bar2, 9))).^(1/2);
fa = Mo .* (2-Mo);
alpha = (3/16)*(5*fa-4);

%% Mach Polynomials
M_l_one_plus = 0.5*(M_l + abs(M_l));
M_r_one_plus = 0.5*(M_r + abs(M_r));

M_l_one_minus = 0.5*(M_l - abs(M_l));
M_r_one_minus = 0.5*(M_r - abs(M_r));

M_l_two_plus = (1/4)*(M_l + 1).^(2);
M_l_two_minus = -(1/4)*(M_l - 1).^(2);

M_r_two_plus = (1/4)*(M_r + 1).^(2);
M_r_two_minus = -(1/4)*(M_r - 1).^(2);

%% Mdot 

%%% Left Plus Mach
M_l_plus = zeros(size(M_l));

one_ind = abs(M_l) >= 1;
M_l_plus(one_ind) = M_l_one_plus(one_ind);
M_l_plus(~one_ind) = M_l_two_plus(~one_ind) .* (1-16.*beta.*M_l_two_minus(~one_ind));

%%% Right Minus Mach
M_r_minus = zeros(size(M_r));

one_ind = abs(M_r) >= 1;
M_r_minus(one_ind) = M_r_one_minus(one_ind);
M_r_minus(~one_ind) = M_r_two_minus(~one_ind) .* (1+16.*beta.*M_r_two_plus(~one_ind));

%%% Mp
Mp = (Kp./fa) .* max(1 - sigma.* M_bar2, 0) .* (P_r - P_l) ./ (rho_half .* c_half.^2);

M_half = M_l_plus + M_r_minus - Mp;
M_wall = M_l_plus + M_r_minus;

%%% Density
rho = rho_r;
rho(M_half > 0) = rho_l(M_half > 0);

mdot_half = c_half .* M_half .* rho;
mdot_wall = c_half .* M_wall .* rho;

%% Pressure
%%% Left Plus
P_l_plus = zeros(size(P_l));

one_ind = abs(M_l) >= 1;
P_l_plus(one_ind) = (1./M_l(one_ind)) .* M_l_one_plus(one_ind);
P_l_plus(~one_ind) = M_l_two_plus(~one_ind) .*...
    ( (2-M_l(~one_ind)) - (16.*alpha(~one_ind).*M_l(~one_ind).*M_l_two_minus(~one_ind)));

%%% Right Plus
P_r_minus = zeros(size(P_r));

one_ind = abs(M_r) >= 1;
P_r_minus(one_ind) = (1./M_r(one_ind)) .* M_r_one_minus(one_ind);
P_r_minus(~one_ind) = M_r_two_minus(~one_ind) .*...
    ( (-2-M_r(~one_ind)) + (16.*alpha(~one_ind).*M_r(~one_ind).*M_l_two_plus(~one_ind)));

%%% Pu
Pu = Ku.*P_l_plus.*P_r_minus.*(rho_l+rho_r).*(fa.*c_half).*(contra_vel_right-contra_vel_left);

P_half = P_l_plus.*P_l + P_r_minus.*P_r - Pu;
P_wall = P_l_plus.*P_l + P_r_minus.*P_r;

%% Flux
mdot = reshape(mdot_half, 1, 1, []);
mdot_bc = reshape(mdot_wall, 1, 1, []);

P = [zeros(1,1,numel(Q_r.q1));...
    reshape(n_x, 1, 1, []) .* reshape(P_half, 1, 1, []);...
    reshape(n_y, 1, 1, []) .* reshape(P_half, 1, 1, []);...
    zeros(1,1,numel(Q_r.q1))];

P_w = [zeros(1,1,numel(Q_r.q1));...
    reshape(n_x, 1, 1, []) .* reshape(P_wall, 1, 1, []);...
    reshape(n_y, 1, 1, []) .* reshape(P_wall, 1, 1, []);...
    zeros(1,1,numel(Q_r.q1))];

f = 0.5.*mdot.*(psi_r+psi_l) - 0.5.*abs(mdot).*(psi_r-psi_l) + P;
f_wall = 0.5.*mdot_bc.*(psi_r+psi_l) + P_w;

f1_wall = reshape(f_wall(1,1,:),size(n_x,1),size(n_x,2));
f2_wall = reshape(f_wall(2,1,:),size(n_x,1),size(n_x,2));
f3_wall = reshape(f_wall(3,1,:),size(n_x,1),size(n_x,2));
f4_wall = reshape(f_wall(4,1,:),size(n_x,1),size(n_x,2));

f1 = reshape(f(1,1,:),size(n_x,1),size(n_x,2));
f2 = reshape(f(2,1,:),size(n_x,1),size(n_x,2));
f3 = reshape(f(3,1,:),size(n_x,1),size(n_x,2));
f4 = reshape(f(4,1,:),size(n_x,1),size(n_x,2));

%% Replace Wall Indices
if strcmp(dir, 'eta')
    % Bottom
    f1(:, 1) = f1_wall(:, 1);
    f2(:, 1) = f2_wall(:, 1);
    f3(:, 1) = f3_wall(:, 1);
    f4(:, 1) = f4_wall(:, 1);

    % Top
    f1(:, end) = f1_wall(:, end);
    f2(:, end) = f2_wall(:, end);
    f3(:, end) = f3_wall(:, end);
    f4(:, end) = f4_wall(:, end);
else
    % Right
    f1(1, :) = f1_wall(1, :);
    f2(1, :) = f2_wall(1, :);
    f3(1, :) = f3_wall(1, :);
    f4(1, :) = f4_wall(1, :);

    % Left
    f1(end, :) = f1_wall(end, :);
    f2(end, :) = f2_wall(end, :);
    f3(end, :) = f3_wall(end, :);
    f4(end, :) = f4_wall(end, :);

end


end