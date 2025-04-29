function [Q_l, Q_r] = muscl(Q, grid, i, j, eps, k, dir, fluid)

arguments
    Q 
    grid 
    i 
    j 
    eps 
    k 
    dir 
    fluid
end

if strcmp(dir, 'xi')
    C1 = 1;
    C2 = 0;

    i_int = i(3:end-2);
    j_int = j(3:end-2);
elseif strcmp(dir, 'eta')
    C1 = 0;
    C2 = 1;

    i_int = i(3:end-2);
    j_int = j(3:end-2);
end

%% Initialize Left/Right States with 1st order
% Boundaries will always be first order to maintain BC's
Q_l.q1 = Q.q1(i-1*C1, j-1*C2);
Q_l.q2 = Q.q2(i-1*C1, j-1*C2);
Q_l.q3 = Q.q3(i-1*C1, j-1*C2);
Q_l.q4 = Q.q4(i-1*C1, j-1*C2);

Q_r.q1 = Q.q1(i, j);
Q_r.q2 = Q.q2(i, j);
Q_r.q3 = Q.q3(i, j);
Q_r.q4 = Q.q4(i, j);

%% Higher Order
fields = ["q1", "q2", "q3", "q4"];
if eps == 1
    [rho_im2, u_im2, v_im2, et_im2, P_im2, T_im2, ht_im2] = Q_to_primitive(Q.q1(i_int - 2*C1, j_int - 2*C2), Q.q2(i_int - 2*C1, j_int - 2*C2),...
        Q.q3(i_int - 2*C1, j_int - 2*C2), Q.q4(i_int - 2*C1, j_int - 2*C2), grid.deltaV(i_int - 2*C1, j_int - 2*C2), fluid);

    [rho_im1, u_im1, v_im1, et_im1, P_im1, T_im1, ht_im1] = Q_to_primitive(Q.q1(i_int - 1*C1, j_int - 1*C2), Q.q2(i_int - 1*C1, j_int - 1*C2),...
        Q.q3(i_int - 1*C1, j_int - 1*C2), Q.q4(i_int - 1*C1, j_int - 1*C2), grid.deltaV(i_int - 1*C1, j_int - 1*C2), fluid);

    [rho_i, u_i, v_i, et_i, P_i, T_i, ht_i] = Q_to_primitive(Q.q1(i_int, j_int), Q.q2(i_int, j_int),...
        Q.q3(i_int, j_int), Q.q4(i_int, j_int), grid.deltaV(i_int, j_int), fluid);

    [rho_ip1, u_ip1, v_ip1, et_ip1, P_ip1, T_ip1, ht_ip1] = Q_to_primitive(Q.q1(i_int + 1*C1, j_int + 1*C2), Q.q2(i_int + 1*C1, j_int + 1*C2),...
        Q.q3(i_int + 1*C1, j_int + 1*C2), Q.q4(i_int + 1*C1, j_int + 1*C2), grid.deltaV(i_int + 1*C1, j_int + 1*C2), fluid);

    [rho_l, rho_r] = interp(rho_im2, rho_im1, rho_i, rho_ip1, k);
    [u_l, u_r] = interp(u_im2, u_im1, u_i, u_ip1, k);
    [v_l, v_r] = interp(v_im2, v_im1, v_i, v_ip1, k);
    [P_l, P_r] = interp(P_im2, P_im1, P_i, P_ip1, k);

    
    Q_l.q1(i_int-1, j_int-1) = rho_l.*grid.deltaV(i_int, j_int);
    Q_r.q1(i_int-1, j_int-1) = rho_r.*grid.deltaV(i_int, j_int);

    Q_l.q2(i_int-1, j_int-1) = rho_l.*u_l.*grid.deltaV(i_int, j_int);
    Q_r.q2(i_int-1, j_int-1) = rho_r.*u_r.*grid.deltaV(i_int, j_int);

    Q_l.q3(i_int-1, j_int-1) = rho_l.*v_l.*grid.deltaV(i_int, j_int);
    Q_r.q3(i_int-1, j_int-1) = rho_r.*v_r.*grid.deltaV(i_int, j_int);

    Q_l.q4(i_int-1, j_int-1) = ( (P_l./(fluid.gamma-1)) + 0.5.*rho_l.*(u_l.^2+v_l.^2) ).*grid.deltaV(i_int, j_int);
    Q_r.q4(i_int-1, j_int-1) = ( (P_r./(fluid.gamma-1)) + 0.5.*rho_r.*(u_r.^2+v_r.^2) ).*grid.deltaV(i_int, j_int);  

end

end    

function [left, right] = interp(im2, im1, i, ip1, k)
im1im2 = im1-im2;
im1im2(im1im2 == 0) = 1e-12;

iim1 = i-im1;
iim1(iim1 == 0) = 1e-12;

ip1i = ip1 - i;
ip1i(ip1i == 0) = 1e-12;

rl = iim1 ./ im1im2;
rr = iim1 ./ ip1i;

phi_l = vanLeer(rl);
phi_r = vanLeer(rr);

phi_l_inv = vanLeer(1./rl);
phi_r_inv = vanLeer(1./rr);

left = im1 +  (1/4) .* ( (1-k)*im1im2.*phi_l + (1+k)*iim1.*phi_l_inv );
right = i - (1/4) .* ( (1+k)*iim1.*phi_r_inv + (1-k)*ip1i.*phi_r );

% left = im1 +  (1/4) .* ( (1-k)*im1im2 + (1+k)*iim1 );
% right = i - (1/4) .* ( (1+k)*iim1 + (1-k)*ip1i );
end

function phi = vanLeer(r)
    phi = (r + abs(r))./(1 + abs(r));
end

function phi = minmod(r)
    phi = max(0, min(1, r));
end