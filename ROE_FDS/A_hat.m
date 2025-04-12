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

%% Eigen Values
lambda1 = (contra_vel-c);
lambda2 = contra_vel;
lambda3 = (contra_vel+c);
lambda4 = contra_vel;

%% Right Eigen Vectors
%%% First Row
r11 = ones(1,numel(u));
r12 = ones(1,numel(u));
r13 = ones(1,numel(u));
r14 = zeros(1,numel(u));

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

%% Left Eigen Vectors
l11 = (2*c.*contra_vel + (gamma-1).*(u.^2+v.^2)) ./ (4*c.^2);
l12 = -(c.*nx + (gamma-1).*u) ./ (2*c.^2);
l13 = (-c.*ny + (gamma-1).*v) ./ (2*c.^2);
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
l44 = zeros(1,numel(u));

%% A hat Matrix

Ahat.a11 = zeros(1,numel(u));
Ahat.a12 = zeros(1,numel(u));
Ahat.a13 = zeros(1,numel(u));
Ahat.a14 = zeros(1,numel(u));
Ahat.a21 = zeros(1,numel(u));
Ahat.a22 = zeros(1,numel(u));
Ahat.a23 = zeros(1,numel(u));
Ahat.a24 = zeros(1,numel(u));
Ahat.a31 = zeros(1,numel(u));
Ahat.a32 = zeros(1,numel(u));
Ahat.a33 = zeros(1,numel(u));
Ahat.a34 = zeros(1,numel(u));
Ahat.a41 = zeros(1,numel(u));
Ahat.a42 = zeros(1,numel(u));
Ahat.a43 = zeros(1,numel(u));
Ahat.a44 = zeros(1,numel(u));
for i = 1:numel(u)
    T_minus = [l11(i) l12(i) l13(i) l14(i);...
               l21(i) l22(i) l23(i) l24(i);...
               l31(i) l32(i) l33(i) l34(i);...
               l41(i) l42(i) l43(i) l44(i)];

    lmb = abs([lambda1(i) 0 0 0;...
               0 lambda2(i) 0 0;...
               0 0 lambda3(i) 0;...
               0 0 0 lambda4(i)]);

    T = [r11(i) r12(i) r13(i) r14(i);...
         r21(i) r22(i) r23(i) r24(i);...
         r31(i) r32(i) r33(i) r34(i);...
         r41(i) r42(i) r43(i) r44(i)];

    A = (1/deltaV(i))*T*lmb*T_minus;

    Ahat.a11(i) = A(1,1);
    Ahat.a12(i) = A(1,2);
    Ahat.a13(i) = A(1,3);
    Ahat.a14(i) = A(1,4);
    Ahat.a21(i) = A(2,1);
    Ahat.a22(i) = A(2,2);
    Ahat.a23(i) = A(2,3);
    Ahat.a24(i) = A(2,4);
    Ahat.a31(i) = A(3,1);
    Ahat.a32(i) = A(3,2);
    Ahat.a33(i) = A(3,3);
    Ahat.a34(i) = A(3,4);
    Ahat.a41(i) = A(4,1);
    Ahat.a42(i) = A(4,2);
    Ahat.a43(i) = A(4,3);
    Ahat.a44(i) = A(4,4);

end


end