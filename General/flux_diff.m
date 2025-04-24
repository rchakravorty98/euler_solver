function [Q_new, E_diff, F_diff] = flux_diff(Q, E, F, grid, deltaT)
Q_new = Q;
nx = grid.nx;
ny = grid.ny;
S_xi = grid.xi.S;
S_eta = grid.eta.S;

%%% E(i+1/2,j) - E(i-1/2,j)
E_diff.e1 = E.e1(2:nx, :).*S_xi(2:nx,:) - E.e1(1:nx-1, :).*S_xi(1:nx-1,:); 
E_diff.e2 = E.e2(2:nx, :).*S_xi(2:nx,:) - E.e2(1:nx-1, :).*S_xi(1:nx-1,:); 
E_diff.e3 = E.e3(2:nx, :).*S_xi(2:nx,:) - E.e3(1:nx-1, :).*S_xi(1:nx-1,:); 
E_diff.e4 = E.e4(2:nx, :).*S_xi(2:nx,:) - E.e4(1:nx-1, :).*S_xi(1:nx-1,:); 

%%% F(i,j+1/2) - F(i,j-1/2)
F_diff.f1 = F.f1(:, 2:ny).*S_eta(:,2:ny) - F.f1(:, 1:ny-1).*S_eta(:, 1:ny-1); 
F_diff.f2 = F.f2(:, 2:ny).*S_eta(:,2:ny) - F.f2(:, 1:ny-1).*S_eta(:, 1:ny-1); 
F_diff.f3 = F.f3(:, 2:ny).*S_eta(:,2:ny) - F.f3(:, 1:ny-1).*S_eta(:, 1:ny-1); 
F_diff.f4 = F.f4(:, 2:ny).*S_eta(:,2:ny) - F.f4(:, 1:ny-1).*S_eta(:, 1:ny-1); 

%%% Q at Next Time Step
Q_new.q1(2:nx, 2:ny) = Q.q1(2:nx, 2:ny) - deltaT.*(E_diff.e1 + F_diff.f1);
Q_new.q2(2:nx, 2:ny) = Q.q2(2:nx, 2:ny) - deltaT.*(E_diff.e2 + F_diff.f2);
Q_new.q3(2:nx, 2:ny) = Q.q3(2:nx, 2:ny) - deltaT.*(E_diff.e3 + F_diff.f3);
Q_new.q4(2:nx, 2:ny) = Q.q4(2:nx, 2:ny) - deltaT.*(E_diff.e4 + F_diff.f4);

end