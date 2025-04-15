function [rho, u, v, et, P, T, ht] = Q_to_primitive(q1, q2, q3, q4, deltaV, fluid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposes Q vector into primitive variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cv = fluid.cv;
R = fluid.R;

rho = q1 ./ deltaV;
u = q2 ./ q1;
v = q3 ./ q1;
et = q4 ./ q1;

T = (et - 0.5*(u.^2+v.^2))/cv;
P = rho .* R .* T; 

ht = et + P./rho;


end