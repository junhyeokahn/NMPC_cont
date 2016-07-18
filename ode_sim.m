function dx_dot = ode_sim(t,x,u,f,B,B_w,w)

dx_dot = f(x) + B*u + B_w*w;

end