function [x_nom, u_nom, ang_accel] = generate_quad_traj(t,Jq,mq,g)

%% Parabola
% om = 2*pi*(1/4);
% x_d = cos(om*t);
% y_d = zeros(length(t),1);
% z_d = -sin(om*t);
% yaw_d = 0.5*cos(om*t);
% % yaw_d = zeros(length(t),1);
% 
% vx_d = -om*sin(om*t);
% vy_d = zeros(length(t),1);
% vz_d = -om*cos(om*t);
% yd_d = -0.5*om*sin(om*t);
% % yd_d = zeros(length(t),1);
% 
% ax_d = -(om^2)*cos(om*t);
% ay_d = zeros(length(t),1);
% az_d =  (om^2)*sin(om*t);
% ydd_d = -0.5*(om^2)*cos(om*t);
% % ydd_d = zeros(length(t),1);
% 
% jx_d =  (om^3)*sin(om*t);
% jy_d = zeros(length(t),1);
% jz_d =  (om^3)*cos(om*t);
% 
% sx_d = (om^4)*cos(om*t);
% sy_d = zeros(length(t),1);
% sz_d =-(om^4)*sin(om*t);

%% Helix
om = 2*pi*0.2;
x_d = sin(om*t);
y_d = cos(om*t);
z_d = -0.1*t;
yaw_d = zeros(length(t),1);

vx_d = om*cos(om*t);
vy_d = -om*sin(om*t);
vz_d = -0.1*ones(length(t),1);
yd_d = zeros(length(t),1);

ax_d = -(om^2)*sin(om*t);
ay_d = -(om^2)*cos(om*t);
az_d = zeros(length(t),1);
ydd_d = zeros(length(t),1);

jx_d = -(om^3)*cos(om*t);
jy_d =  (om^3)*sin(om*t);
jz_d = zeros(length(t),1);

sx_d = (om^4)*sin(om*t);
sy_d = (om^4)*cos(om*t);
sz_d = zeros(length(t),1);

%% Plot
figure()
plot3(x_d,y_d,z_d); grid on
xlabel('x'); ylabel('y'); zlabel('h');
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');


%% Mapping

%Get control & nominal state
u_nom = zeros(length(t),4);
x_nom = zeros(length(t),14);
ang_accel = zeros(length(t),3);

x_nom(:,1:6) = [x_d,y_d,z_d,vx_d,vy_d,vz_d];
x_nom(:,11) = yaw_d;

for i = 1:length(t)
   th_vect = -[ax_d(i);ay_d(i);az_d(i)-g];
   zb_ax = th_vect/norm(th_vect);
   
   thrust = mq*norm(th_vect);
   thrust_d = -mq*(zb_ax'*[jx_d(i);jy_d(i);jz_d(i)]);
   x_nom(i,7) = thrust; %thrust
   x_nom(i,8) = thrust_d ; %thrust_dot
   
   beta_a = -cos(yaw_d(i))*ax_d(i) - sin(yaw_d(i))*ay_d(i);
   beta_b = g- az_d(i);
   beta_c = -sin(yaw_d(i))*ax_d(i) + cos(yaw_d(i))*ay_d(i);
      
   roll = atan2(beta_c,norm([beta_a;beta_b]));%roll
   pitch = atan2(beta_a,beta_b);%pitch
      
   x_nom(i,9) = roll; 
   x_nom(i,10) = pitch;
   
   xb_ax = [cos(pitch)*cos(yaw_d(i));
            cos(pitch)*sin(yaw_d(i));
            -sin(pitch)];
   yb_ax = cross(zb_ax,xb_ax);
   
   h_om =-(1/thrust)*(mq*[jx_d(i);jy_d(i);jz_d(i)] + ...
                      thrust_d*zb_ax);
           
   x_nom(i,12) = -h_om'*yb_ax; %p
   x_nom(i,13) = h_om'*xb_ax; %q
   x_nom(i,14) = -x_nom(i,13)*tan(roll) + ...
              yd_d(i)*cos(pitch)*(tan(roll)*sin(roll) + cos(roll)); %r
   
   thrust_dd = mq*(-zb_ax'*[sx_d(i);sy_d(i);sz_d(i)]) + thrust*(h_om'*h_om); %thrust_ddot      
   h_al =-(1/thrust)*(mq*[sx_d(i);sy_d(i);sz_d(i)] + thrust_dd*zb_ax + 2*thrust_d*h_om + ...
                      thrust*cross(x_nom(i,12:14)',h_om));

   pd = -h_al'*yb_ax;
   qd = h_al'*xb_ax;
   
   phi_d = x_nom(i,12) + yd_d(i)*sin(pitch);
   theta_d = sec(roll)*(x_nom(i,13) - yd_d(i)*sin(roll)*cos(pitch));
   
   rd = -qd*tan(roll) - x_nom(i,13)*phi_d*(sec(roll))^2 + ...
         (ydd_d(i)*cos(pitch) - yd_d(i)*theta_d*sin(pitch))*(tan(roll)*sin(roll) + cos(roll)) + ...
         yd_d(i)*phi_d*cos(pitch)*(sin(roll)*(sec(roll))^2 + tan(roll)*cos(roll) - sin(roll));
   
   ang_accel(i,:) = [pd,qd,rd];
   torque = Jq*[pd;qd;rd] + Skew(x_nom(i,12:14)')*Jq*x_nom(i,12:14)';
   
   u_nom(i,:) = [thrust_dd, torque'];      
    
end

end