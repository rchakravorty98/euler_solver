function Ahat = A_hat(u,v,c,gamma,Sx,Sy,S, deltaV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the A hat matrix in Roe FDS formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    u       % u-component [m/s]
    v       % v-Component [m/s]
    c       % Speed of Sound [m/s]
    gamma   % Specific Heat Ratio
    Sx      % Projected Cell Face Area in X
    Sy      % Projected Cell Face Area in y
    S       % Cell Face Area
    deltaV  % Cell Volume
end

nx = Sx./S;
ny = Sy./S;
contra_vel = u.*nx + v.*ny;     % Contravariant Velocity
h = c.^2./(gamma-1);            % Enthalphy
ek = 0.5.*(u.^2+v.^2);          % Kinetic Energy

num_elem = numel(nx);

%% Eigen Values
lambda1 = (contra_vel-c);
lambda2 = contra_vel;
lambda3 = (contra_vel+c);
lambda4 = contra_vel;

lmb_array = abs([reshape(lambda1, 1, 1,[]) zeros(1,1,num_elem) zeros(1,1,num_elem) zeros(1,1,num_elem);...
    zeros(1,1,num_elem) reshape(lambda2, 1, 1,[])  zeros(1,1,num_elem) zeros(1,1,num_elem);...
    zeros(1,1,num_elem) zeros(1,1,num_elem) reshape(lambda3, 1, 1,[])  zeros(1,1,num_elem);...
    zeros(1,1,num_elem) zeros(1,1,num_elem) zeros(1,1,num_elem) reshape(lambda4, 1, 1,[])]);

%% Right Eigen Vectors
%%% First Row
r11 = ones(size(nx));
r12 = ones(size(nx));
r13 = ones(size(nx));
r14 = zeros(size(nx));

%%% Second Row
r21 = u-c.*nx;
r22 = u;
r23 = u+c.*nx;
r24 = ny;

%%% Third Row
r31 = v-c.*ny;
r32 = v;
r33 = v+c.*ny;
r34 = -nx;

%%% Fourth Row
r41 = h + ek - c.* contra_vel;
r42 = ek;
r43 = h + ek + c.* contra_vel;
r44 = u.*ny - v.*nx;

right_vec = [reshape(r11,1,1,[]) reshape(r12,1,1,[]) reshape(r13,1,1,[]) reshape(r14,1,1,[]);...
    reshape(r21,1,1,[]) reshape(r22,1,1,[]) reshape(r23,1,1,[]) reshape(r24,1,1,[]);...
    reshape(r31,1,1,[]) reshape(r32,1,1,[]) reshape(r33,1,1,[]) reshape(r34,1,1,[]);...
    reshape(r41,1,1,[]) reshape(r42,1,1,[]) reshape(r43,1,1,[]) reshape(r44,1,1,[])];

%% Left Eigen Vectors
l11 = (2*c.*contra_vel + (gamma-1).*(u.^2+v.^2)) ./ (4*c.^2);
l12 = -(c.*nx + (gamma-1).*u) ./ (2*c.^2);
l13 = -(c.*ny + (gamma-1).*v) ./ (2*c.^2);
l14 = (gamma-1) ./ (2*c.^2);

l21 = 1 - ((gamma-1).*(u.^2+v.^2)) ./ (2*c.^2);
l22 = ((gamma-1)*u) ./ (c.^2);
l23 = ((gamma-1)*v) ./ (c.^2);
l24 = -(gamma-1) ./ (c.^2);

l31 = (-2*c.*contra_vel + (gamma-1).*(u.^2+v.^2)) ./ (4*c.^2);
l32 = (c.*nx - (gamma-1).*u) ./ (2*c.^2);
l33 = (c.*ny - (gamma-1).*v) ./ (2*c.^2);
l34 = (gamma-1) ./ (2*c.^2);

l41 = -u.*ny + v.*nx;
l42 = ny;
l43 = -nx;
l44 = zeros(size(nx));

left_vec = [reshape(l11,1,1,[]) reshape(l12,1,1,[]) reshape(l13,1,1,[]) reshape(l14,1,1,[]);...
    reshape(l21,1,1,[]) reshape(l22,1,1,[]) reshape(l23,1,1,[]) reshape(l24,1,1,[]);...
    reshape(l31,1,1,[]) reshape(l32,1,1,[]) reshape(l33,1,1,[]) reshape(l34,1,1,[]);...
    reshape(l41,1,1,[]) reshape(l42,1,1,[]) reshape(l43,1,1,[]) reshape(l44,1,1,[])];

%% A hat Matrix

temp = pagemtimes(right_vec, lmb_array);
Ahat_temp = pagemtimes(temp, left_vec);

Ahat.a11 = reshape(Ahat_temp(1,1,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a12 = reshape(Ahat_temp(1,2,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a13 = reshape(Ahat_temp(1,3,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a14 = reshape(Ahat_temp(1,4,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a21 = reshape(Ahat_temp(2,1,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a22 = reshape(Ahat_temp(2,2,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a23 = reshape(Ahat_temp(2,3,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a24 = reshape(Ahat_temp(2,4,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a31 = reshape(Ahat_temp(3,1,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a32 = reshape(Ahat_temp(3,2,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a33 = reshape(Ahat_temp(3,3,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a34 = reshape(Ahat_temp(3,4,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a41 = reshape(Ahat_temp(4,1,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a42 = reshape(Ahat_temp(4,2,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a43 = reshape(Ahat_temp(4,3,:),size(u,1),size(u,2)) ./ deltaV;
Ahat.a44 = reshape(Ahat_temp(4,4,:),size(u,1),size(u,2)) ./ deltaV;

end