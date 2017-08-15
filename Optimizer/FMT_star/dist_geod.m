function dist = dist_geod(q_start,q_end,EPS)

T = norm(q_end.coord-q_start.coord)/1.7;
dist = q_start.cost*exp(-EPS.lambda*T) + EPS.d_bar*(1-exp(-EPS.lambda*T));

end