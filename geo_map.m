function X_dot = geo_map(x_start,x_end,xi_start,xi_end,phi_d,n)

xi_dot = xi_end - xi_start;

X_dot = zeros(n,2);

X_dot(:,1) = phi_d(x_start)\xi_dot;
X_dot(:,2) = phi_d(x_end)\xi_dot;

end