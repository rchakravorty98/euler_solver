function [E, F] = roe_fds(Q, grid, fluid)

nx = grid.nx;
ny = grid.ny;

R = fluid.R;            % Ideal Gas Constant [J/kg*K]
cp = fluid.cp;          % Specific Heat [J/kg*k]
gamma = fluid.gamma;    % Specific Heat Ratio
cv = fluid.cv; 

%% xi direction
j_idx = 2:ny;   % Need fluxes over primary grid

E.e1 = zeros(numel(2:nx+1), numel(2:ny));
E.e2 = zeros(numel(2:nx+1), numel(2:ny));
E.e3 = zeros(numel(2:nx+1), numel(2:ny));
E.e4 = zeros(numel(2:nx+1), numel(2:ny));
for i = 2:nx+1
    %%% Left State
    [rho_l, u_l, v_l, et_l, P_l, T_l] = Q_to_primitive(Q.q1(i-1,j_idx), Q.q2(i-1,j_idx),...
        Q.q3(i-1,j_idx), Q.q4(i-1,j_idx), grid.deltaV(i-1,j_idx), fluid);

    ht_l = et_l + P_l./rho_l;                           % Left Total Enthalpy

    %%% Right State
    [rho_r, u_r, v_r, et_r, P_r, T_r] = Q_to_primitive(Q.q1(i,j_idx), Q.q2(i,j_idx),...
        Q.q3(i,j_idx), Q.q4(i,j_idx), grid.deltaV(i,j_idx), fluid);

    ht_r = et_r + P_r./rho_r;                       % Right Total Enthalpy

    %%% Roe Averages
    u_bar = ((rho_r.^(1/2) .* u_r) + (rho_l.^(1/2) .* u_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    v_bar = ((rho_r.^(1/2) .* v_r) + (rho_l.^(1/2) .* v_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    rho_bar = (rho_r .* rho_l).^(1/2);
    ht_bar = ((rho_r.^(1/2) .* ht_r) + (rho_l.^(1/2) .* ht_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    c_bar = ((gamma-1) * (ht_bar - 0.5*(u_bar.^2 + v_bar.^2))).^(1/2);
    
    %%% Jacobian
    Ahat = A_hat(u_bar, v_bar, c_bar,...
        gamma, grid.xi.Sx(i-1,:), grid.xi.Sy(i-1,:),...
        grid.xi.S(i-1,:), grid.deltaV(i,j_idx));

    %%% Construct Flux
    for j = 1:numel(j_idx)
        A_hat_ij = [Ahat.a11(j) Ahat.a12(j) Ahat.a13(j) Ahat.a14(j);
            Ahat.a21(j)  Ahat.a22(j) Ahat.a23(j) Ahat.a24(j);
            Ahat.a31(j)  Ahat.a32(j) Ahat.a33(j) Ahat.a34(j);
            Ahat.a41(j)  Ahat.a42(j) Ahat.a43(j) Ahat.a44(j)];

        n_x = grid.xi.Sx(i-1,j)/grid.xi.S(i-1,j);
        n_y = grid.xi.Sy(i-1,j)/grid.xi.S(i-1,j);

        %%% Left
        contra_vel_left = u_l(j)* n_x + v_l(j)*n_y;

        E_l = [rho_l(j)*contra_vel_left;...
            rho_l(j)*u_l(j)*contra_vel_left + n_x * P_l(j);...
            rho_l(j)*v_l(j)*contra_vel_left + n_y * P_l(j);...
            rho_l(j)*ht_l(j)*contra_vel_left];

        Q_l = [Q.q1(i-1,j+1);...
            Q.q2(i-1,j+1);...
            Q.q3(i-1,j+1);...
            Q.q4(i-1,j+1)];

        %%% Right
        contra_vel_right = u_r(j)* n_x + v_r(j)*n_y;

        E_r = [rho_r(j)*contra_vel_right;...
            rho_r(j)*u_r(j)*contra_vel_right + n_x * P_r(j);...
            rho_r(j)*v_r(j)*contra_vel_right + n_y * P_r(j);...
            rho_r(j)*ht_r(j)*contra_vel_right];
        
        Q_r = [Q.q1(i,j+1);...
            Q.q2(i,j+1);...
            Q.q3(i,j+1);...
            Q.q4(i,j+1)];

        E_half = 0.5 * (E_l + E_r) - 0.5*A_hat_ij*(Q_r - Q_l);

        E.e1(i-1,j) = E_half(1);
        E.e2(i-1,j) = E_half(2);
        E.e3(i-1,j) = E_half(3);
        E.e4(i-1,j) = E_half(4);


    end

end

%% eta direction
i_idx = 2:nx;   % Need fluxes over primary grid

F.f1 = zeros(numel(2:nx), numel(2:ny+1));
F.f2 = zeros(numel(2:nx), numel(2:ny+1));
F.f3 = zeros(numel(2:nx), numel(2:ny+1));
F.f4 = zeros(numel(2:nx), numel(2:ny+1));
for j = 2:ny+1
    %%% Left State
    [rho_l, u_l, v_l, et_l, P_l, T_l] = Q_to_primitive(Q.q1(i_idx,j-1), Q.q2(i_idx,j-1),...
        Q.q3(i_idx,j-1), Q.q4(i_idx,j-1), grid.deltaV(i_idx,j-1), fluid);

    ht_l = et_l + P_l./rho_l;                           % Left Total Enthalpy

    %%% Right State
    [rho_r, u_r, v_r, et_r, P_r, T_r] = Q_to_primitive(Q.q1(i_idx,j), Q.q2(i_idx,j),...
        Q.q3(i_idx,j), Q.q4(i_idx,j), grid.deltaV(i_idx,j), fluid);

    ht_r = et_r + P_r./rho_r;                       % Right Total Enthalpy

    %%% Roe Averages
    u_bar = ((rho_r.^(1/2) .* u_r) + (rho_l.^(1/2) .* u_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    v_bar = ((rho_r.^(1/2) .* v_r) + (rho_l.^(1/2) .* v_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    rho_bar = (rho_r .* rho_l).^(1/2);
    ht_bar = ((rho_r.^(1/2) .* ht_r) + (rho_l.^(1/2) .* ht_l)) ./ ...
        (rho_r.^(1/2) + rho_l.^(1/2));
    c_bar = ((gamma-1) * (ht_bar - 0.5*(u_bar.^2 + v_bar.^2))).^(1/2);
    
    %%% Jacobian
    Ahat = A_hat(u_bar, v_bar, c_bar,...
        gamma, grid.eta.Sx(i_idx-1,:), grid.eta.Sy(i_idx-1,:),...
        grid.eta.S(i_idx-1,:), grid.deltaV(i_idx,j));

    
    %%% Construct Flux
    for i = 1:numel(i_idx)
        A_hat_ij = [Ahat.a11(i) Ahat.a12(i) Ahat.a13(i) Ahat.a14(i);
            Ahat.a21(i)  Ahat.a22(i) Ahat.a23(i) Ahat.a24(i);
            Ahat.a31(i)  Ahat.a32(i) Ahat.a33(i) Ahat.a34(i);
            Ahat.a41(i)  Ahat.a42(i) Ahat.a43(i) Ahat.a44(i)];

        n_x = grid.eta.Sx(i,j-1)/grid.eta.S(i,j-1);
        n_y = grid.eta.Sy(i,j-1)/grid.eta.S(i,j-1);

        %%% Left
        contra_vel_left = u_l(i)* n_x + v_l(i)*n_y;

        F_l = [rho_l(i)*contra_vel_left;...
            rho_l(i)*u_l(i)*contra_vel_left + n_x * P_l(i);...
            rho_l(i)*v_l(i)*contra_vel_left + n_y * P_l(i);...
            rho_l(i)*ht_l(i)*contra_vel_left];

        Q_l = [Q.q1(i+1,j-1);...
            Q.q2(i+1,j-1);...
            Q.q3(i+1,j-1);...
            Q.q4(i+1,j-1)];

        %%% Right
        contra_vel_right = u_r(j)* n_x + v_r(j)*n_y;

        F_r = [rho_r(i)*contra_vel_right;...
            rho_r(i)*u_r(i)*contra_vel_right + n_x * P_r(i);...
            rho_r(i)*v_r(i)*contra_vel_right + n_y * P_r(i);...
            rho_r(i)*ht_r(i)*contra_vel_right];
        
        Q_r = [Q.q1(i+1,j);...
            Q.q2(i+1,j);...
            Q.q3(i+1,j);...
            Q.q4(i+1,j)];

        F_half = 0.5 * (F_l + F_r) - 0.5*A_hat_ij*(Q_r - Q_l);

        F.f1(i,j-1) = F_half(1);
        F.f2(i,j-1) = F_half(2);
        F.f3(i,j-1) = F_half(3);
        F.f4(i,j-1) = F_half(4);


    end

end



end