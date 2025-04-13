function res = residual(Q, Q_new, grid, free_stream)

q1_diff = abs(Q_new.q1 - Q.q1) ./ grid.deltaV;
q2_diff = abs((Q_new.q2./Q_new.q1) - (Q.q2./Q.q1));
q3_diff = abs((Q_new.q3./Q_new.q1) - (Q.q3./Q.q1));
q4_diff = abs((Q_new.q4./Q_new.q1) - (Q.q4./Q.q1));

q1_norm = q1_diff/free_stream.rho_ref;
q2_norm = q2_diff/free_stream.u_ref;
q3_norm = q3_diff/free_stream.u_ref;
q4_norm = q4_diff/(free_stream.rho_et/free_stream.rho_ref);

res = max(max([q1_norm q2_norm q3_norm q4_norm]));


end