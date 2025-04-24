function [Q_l, Q_r] = muscl(Q, grid, i, j, eps, k, dir)

arguments
    Q 
    grid 
    i 
    j 
    eps 
    k 
    dir 
end

if strcmp(dir, 'xi')
    C1 = 1;
    C2 = 0;

    i_int = i(2:end-1);
    j_int = j(2:end-1);
elseif strcmp(dir, 'eta')
    C1 = 0;
    C2 = 1;

    i_int = i(2:end-1);
    j_int = j(2:end-1);
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
    for field = fields
        q = Q.(field);

        q_im2 = q(i_int - 2*C1, j_int - 2*C2);
        q_im1 = q(i_int - 1*C1, j_int - 1*C2);
        q_i   = q(i_int,        j_int);
        q_ip1 = q(i_int + 1*C1, j_int + 1*C2);

        dq_im1_im2 = q_im1-q_im2;
        dq_i_im1 = q_i-q_im1;
        dq_ip1_i = q_ip1 - q_i;
        
        rl = (q_i-q_im1)./(q_im1-q_im2);
        rr = (q_i-q_im1)./(q_ip1-q_i);

        phi_l = vanLeer(rl);
        phi_r = vanLeer(rr);

        phi_l_inv = vanLeer(1./rl);
        phi_r_inv = vanLeer(1./rr);

        Q_l.(field)(i_int-1, j_int-1) = q_im1 +  (1/4) .* ( (1-k)*dq_im1_im2.*phi_l + (1+k)*dq_i_im1.*phi_l_inv );
        Q_r.(field)(i_int-1, j_int-1) = q_i - (1/4) .* ( (1+k)*dq_i_im1.*phi_r_inv + (1-k)*dq_ip1_i.*phi_r );

        % Q_l.(field)(i_int-1, j_int-1) = q_im1 + (1/4) .* ( (1-k)*dq_im1_im2 + (1+k)*dq_i_im1 );
        % Q_r.(field)(i_int-1, j_int-1) = q_i - (1/4) .* ( (1+k)*dq_i_im1 + (1-k)*dq_ip1_i );

    end
end

end    

function phi = vanLeer(r)
    phi = (r + abs(r))./(1 + abs(r) + 1e-12);
    phi(isnan(phi)) = 0;
end

function phi = minmod(r)
    phi = max(0, min(1,r));
    phi(isnan(phi)) = 0;
end

function phi = MC(r)
    C1 = min(2*r, 0.5*(1+r));
    C2 = min(C1, 2);
    phi = max(0, C2);
end