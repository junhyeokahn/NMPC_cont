function dx_dot = quad_ode(t,x,t_span,u_nom,k,f,B,B_w,w)

global kp_om;

uc = interp1(t_span,u_nom,t) + k; %thrust_dot,roll_dot,pitch_dot,yaw_dot

euler_dot_des = uc(2:4)';
om_des = R_om(x(7:9))*euler_dot_des;

M = kp_om*(om_des - x(10:12));

u = [uc(1);M];

dx_dot = f(x) + B*u + B_w*w;

end