function [E, F] = roe_fds(Q, grid, fluid)
eps = 0;
k = 0.5;

nx = grid.nx;
ny = grid.ny;

R = fluid.R;            % Ideal Gas Constant [J/kg*K]
cp = fluid.cp;          % Specific Heat [J/kg*k]
gamma = fluid.gamma;    % Specific Heat Ratio
cv = fluid.cv; 

%% xi direction
i_idx = 2:nx+1;
j_idx = 2:ny;   % Need fluxes over primary grid

E.e1 = zeros(numel(2:nx+1), numel(2:ny));
E.e2 = zeros(numel(2:nx+1), numel(2:ny));
E.e3 = zeros(numel(2:nx+1), numel(2:ny));
E.e4 = zeros(numel(2:nx+1), numel(2:ny));

[Q_l, Q_r] = muscl(Q, grid, i_idx, j_idx, eps, k, 'xi');

%%% Left State
[rho_l, u_l, v_l, et_l, P_l, T_l, ht_l] = Q_to_primitive(Q_l.q1, Q_l.q2,...
    Q_l.q3, Q_l.q4, grid.deltaV(i_idx-1,j_idx), fluid);

%%% Right State
[rho_r, u_r, v_r, et_r, P_r, T_r, ht_r] = Q_to_primitive(Q_r.q1, Q_r.q2,...
    Q_r.q3, Q_r.q4, grid.deltaV(i_idx,j_idx), fluid);

%%% Roe Averages
[u_bar, v_bar, ht_bar, c_bar] = roe_avg(rho_l, rho_r, u_l, u_r, v_l, v_r, ht_l, ht_r, fluid);

%%% Jacobian
deltaV_avg = (grid.deltaV(2:nx+1,j_idx)+ grid.deltaV(1:nx,j_idx))/2;
Ahat = A_hat(u_bar, v_bar, c_bar,...
    gamma, grid.xi.Sx, grid.xi.Sy,...
    grid.xi.S, deltaV_avg);

Ahat_re = [reshape(Ahat.a11,1,1,[]) reshape(Ahat.a12,1,1,[]) reshape(Ahat.a13,1,1,[]) reshape(Ahat.a14,1,1,[]);...
    reshape(Ahat.a21,1,1,[]) reshape(Ahat.a22,1,1,[]) reshape(Ahat.a23,1,1,[]) reshape(Ahat.a24,1,1,[]);...
    reshape(Ahat.a31,1,1,[]) reshape(Ahat.a32,1,1,[]) reshape(Ahat.a33,1,1,[]) reshape(Ahat.a34,1,1,[]);...
    reshape(Ahat.a41,1,1,[]) reshape(Ahat.a42,1,1,[]) reshape(Ahat.a43,1,1,[]) reshape(Ahat.a44,1,1,[])];

n_x = grid.xi.Sx./grid.xi.S;
n_y = grid.xi.Sy./grid.xi.S;

%%% Left
contra_vel_left = u_l.* n_x + v_l.*n_y;

E_l = [reshape(rho_l,1,1,[]) .* reshape(contra_vel_left,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(u_l,1,1,[]).*reshape(contra_vel_left,1,1,[]) + reshape(n_x,1,1,[]) .* reshape(P_l,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(v_l,1,1,[]).*reshape(contra_vel_left,1,1,[]) + reshape(n_y,1,1,[]) .* reshape(P_l,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(ht_l,1,1,[]).*reshape(contra_vel_left,1,1,[])];

Q_l = [reshape(Q_l.q1,1,1,[]);...
       reshape(Q_l.q2,1,1,[]);...
       reshape(Q_l.q3,1,1,[]);...
       reshape(Q_l.q4,1,1,[])];

%%% Right
contra_vel_right = u_r.* n_x + v_r.*n_y;

E_r = [reshape(rho_r,1,1,[]) .* reshape(contra_vel_right,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(u_r,1,1,[]).*reshape(contra_vel_right,1,1,[]) + reshape(n_x,1,1,[]) .* reshape(P_r,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(v_r,1,1,[]).*reshape(contra_vel_right,1,1,[]) + reshape(n_y,1,1,[]) .* reshape(P_r,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(ht_r,1,1,[]).*reshape(contra_vel_right,1,1,[])];

Q_r = [reshape(Q_r.q1,1,1,[]);...
       reshape(Q_r.q2,1,1,[]);...
       reshape(Q_r.q3,1,1,[]);...
       reshape(Q_r.q4,1,1,[])];

E_half = 0.5 * (E_l + E_r) - 0.5*pagemtimes(Ahat_re,(Q_r - Q_l));

E.e1 = reshape(E_half(1,1,:),size(n_x,1),size(n_x,2));
E.e2 = reshape(E_half(2,1,:),size(n_x,1),size(n_x,2));
E.e3 = reshape(E_half(3,1,:),size(n_x,1),size(n_x,2));
E.e4 = reshape(E_half(4,1,:),size(n_x,1),size(n_x,2));

%% eta direction
i_idx = 2:nx;   % Need fluxes over primary grid
j_idx = 2:ny+1;

F.f1 = zeros(numel(2:nx), numel(2:ny+1));
F.f2 = zeros(numel(2:nx), numel(2:ny+1));
F.f3 = zeros(numel(2:nx), numel(2:ny+1));
F.f4 = zeros(numel(2:nx), numel(2:ny+1));

[Q_l, Q_r] = muscl(Q, grid, i_idx, j_idx, eps, k, 'eta');

%%% Left State
[rho_l, u_l, v_l, et_l, P_l, T_l, ht_l] = Q_to_primitive(Q_l.q1, Q_l.q2,...
    Q_l.q3, Q_l.q4, grid.deltaV(i_idx,j_idx-1), fluid);

%%% Right State
[rho_r, u_r, v_r, et_r, P_r, T_r, ht_r] = Q_to_primitive(Q_r.q1, Q_r.q2,...
    Q_r.q3, Q_r.q4, grid.deltaV(i_idx,j_idx), fluid);

%%% Roe Averages
[u_bar, v_bar, ht_bar, c_bar] = roe_avg(rho_l, rho_r, u_l, u_r, v_l, v_r, ht_l, ht_r, fluid);

%%% Jacobian
deltaV_avg = (grid.deltaV(i_idx,2:ny+1)+ grid.deltaV(i_idx,1:ny))/2;
Ahat = A_hat(u_bar, v_bar, c_bar,...
    gamma, grid.eta.Sx, grid.eta.Sy,...
    grid.eta.S, deltaV_avg);

Ahat_re = [reshape(Ahat.a11,1,1,[]) reshape(Ahat.a12,1,1,[]) reshape(Ahat.a13,1,1,[]) reshape(Ahat.a14,1,1,[]);...
    reshape(Ahat.a21,1,1,[]) reshape(Ahat.a22,1,1,[]) reshape(Ahat.a23,1,1,[]) reshape(Ahat.a24,1,1,[]);...
    reshape(Ahat.a31,1,1,[]) reshape(Ahat.a32,1,1,[]) reshape(Ahat.a33,1,1,[]) reshape(Ahat.a34,1,1,[]);...
    reshape(Ahat.a41,1,1,[]) reshape(Ahat.a42,1,1,[]) reshape(Ahat.a43,1,1,[]) reshape(Ahat.a44,1,1,[])];

n_x = grid.eta.Sx./grid.eta.S;
n_y = grid.eta.Sy./grid.eta.S;

%%% Left
contra_vel_left = u_l.* n_x + v_l.*n_y;

F_l = [reshape(rho_l,1,1,[]) .* reshape(contra_vel_left,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(u_l,1,1,[]).*reshape(contra_vel_left,1,1,[]) + reshape(n_x,1,1,[]) .* reshape(P_l,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(v_l,1,1,[]).*reshape(contra_vel_left,1,1,[]) + reshape(n_y,1,1,[]) .* reshape(P_l,1,1,[]);...
    reshape(rho_l,1,1,[]).*reshape(ht_l,1,1,[]).*reshape(contra_vel_left,1,1,[])];

Q_l = [reshape(Q_l.q1,1,1,[]);...
       reshape(Q_l.q2,1,1,[]);...
       reshape(Q_l.q3,1,1,[]);...
       reshape(Q_l.q4,1,1,[])];

%%% Right
contra_vel_right = u_r.* n_x + v_r.*n_y;

F_r = [reshape(rho_r,1,1,[]) .* reshape(contra_vel_right,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(u_r,1,1,[]).*reshape(contra_vel_right,1,1,[]) + reshape(n_x,1,1,[]) .* reshape(P_r,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(v_r,1,1,[]).*reshape(contra_vel_right,1,1,[]) + reshape(n_y,1,1,[]) .* reshape(P_r,1,1,[]);...
    reshape(rho_r,1,1,[]).*reshape(ht_r,1,1,[]).*reshape(contra_vel_right,1,1,[])];

Q_r = [reshape(Q_r.q1,1,1,[]);...
       reshape(Q_r.q2,1,1,[]);...
       reshape(Q_r.q3,1,1,[]);...
       reshape(Q_r.q4,1,1,[])];

F_half = 0.5 * (F_l + F_r) - 0.5*pagemtimes(Ahat_re,(Q_r - Q_l));

F.f1 = reshape(F_half(1,1,:),size(n_x,1),size(n_x,2));
F.f2 = reshape(F_half(2,1,:),size(n_x,1),size(n_x,2));
F.f3 = reshape(F_half(3,1,:),size(n_x,1),size(n_x,2));
F.f4 = reshape(F_half(4,1,:),size(n_x,1),size(n_x,2));

end