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
elseif strcmp(dir, 'eta')
    C1 = 0;
    C2 = 1;
end


nx = grid.nx;
ny = grid.ny;

if eps == 0
    Q_l.q1 = Q.q1(i-1*C1, j-1*C2);
    Q_l.q2 = Q.q2(i-1*C1, j-1*C2);
    Q_l.q3 = Q.q3(i-1*C1, j-1*C2);
    Q_l.q4 = Q.q4(i-1*C1, j-1*C2);

    Q_r.q1 = Q.q1(i, j);
    Q_r.q2 = Q.q2(i, j);
    Q_r.q3 = Q.q3(i, j);
    Q_r.q4 = Q.q4(i, j);

end



end