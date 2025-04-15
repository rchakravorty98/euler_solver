function [u_bar, v_bar, ht_bar, c_bar] = roe_avg(rho_l, rho_r, u_l, u_r, v_l, v_r, ht_l, ht_r, fluid)
gamma = fluid.gamma;

%%% Velocity
u_bar = ((rho_r.^(1/2) .* u_r) + (rho_l.^(1/2) .* u_l)) ./ ...
    (rho_r.^(1/2) + rho_l.^(1/2));
v_bar = ((rho_r.^(1/2) .* v_r) + (rho_l.^(1/2) .* v_l)) ./ ...
    (rho_r.^(1/2) + rho_l.^(1/2));

%%% Enthalpy
ht_bar = ((rho_r.^(1/2) .* ht_r) + (rho_l.^(1/2) .* ht_l)) ./ ...
    (rho_r.^(1/2) + rho_l.^(1/2));

%%% Speed of Sound
c_bar = ((gamma-1) * (ht_bar - 0.5*(u_bar.^2 + v_bar.^2))).^(1/2);

end