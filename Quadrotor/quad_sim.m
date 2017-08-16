function dxth_dot = quad_sim(t,x_th,t_span,u_dyn_nom,u_dyn_fb,f,B,B_w,w)

global Jq kp_om ki_om;

%Actual states
x = x_th(1:12);
q = x(7:9);
om = x(10:12);
eul_rate = R_eul(q)*om;

%Internal states
th = x_th(13);
rate_err_int = x_th(14:16);

%Dynamic control
%Interpolate for nominal, keep feedback in ZOH
u_dyn = interp1(t_span,u_dyn_nom,t) + u_dyn_fb; %[th_dot, eul_rate]
om_des = R_om(q)*u_dyn(2:4)';

%PI controller for rate
M = cross(om,Jq*om) + kp_om*(om_des-om) + ki_om*R_om(q)*rate_err_int;

%Actual control
u = [th; M];

%Integrate dynamics
dx_dot = f(x) + B(x)*u + B_w*w;

%Integrate dynamic controller 'internal states'
dth_dot = u_dyn(1);
drate_err_int_dot = u_dyn(2:4)'-eul_rate;

dxth_dot = [dx_dot; dth_dot; drate_err_int_dot];

end