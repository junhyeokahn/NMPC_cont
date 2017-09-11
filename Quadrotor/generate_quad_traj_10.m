function [t_vec,x_nom, u_nom, th_nom, accel] = generate_quad_traj_10(dt,Jq,mq,g,poly_file)

sim_case = 'hover';

switch (sim_case)
    case 'fig8'       
        T = 10;
        t_vec = (0:dt:T)';
        om = 2*pi*0.125;
        x_d = 2*cos(om*t_vec);
        y_d = 3*0.5*sin(2*om*tt_vec);
        z_d = zeros(length(t_vec),1);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = -2*om*sin(om*t_vec);
        vy_d =  3*om*cos(2*om*t_vec);
        vz_d = zeros(length(t_vec),1);
        yd_d = zeros(length(t_vec,1));
        
        ax_d = -2*(om^2)*cos(om*t_vec);
        ay_d = -3*2*(om^2)*sin(2*om*t_vec);
        az_d =  zeros(length(t_vec),1);
        
        jx_d =  2*(om^3)*sin(om*t_vec);
        jy_d = -3*4*(om^3)*cos(2*om*t_vec);
        jz_d =  zeros(length(t_vec),1);
        
    case 'circle'
        T = 20;
        t_vec = (0:dt:T)';
        radius = 1;
        
        om = 2*pi*(1.5/10);
        x_d = radius*cos(om*t_vec);
%         z_d = -sin(om*t_vec);
        z_d = -0.5*ones(length(t_vec),1);
        y_d = radius*sin(om*t_vec);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = -radius*om*sin(om*t_vec);
%         vz_d = -om*cos(om*t_vec);
        vz_d = zeros(length(t_vec),1);
        vy_d = radius*om*cos(om*t_vec);
        yd_d = zeros(length(t_vec),1);
        
        ax_d = -radius*(om^2)*cos(om*t_vec);
%         az_d = (om^2)*sin(om*t_vec);
        az_d = zeros(length(t_vec),1);
        ay_d = -radius*(om^2)*sin(om*t_vec);
        
        jx_d =  radius*(om^3)*sin(om*t_vec);
%         jz_d = (om^3)*cos(om*t_vec);
        jz_d = zeros(length(t_vec),1);
        jy_d = -radius*(om^3)*cos(om*t_vec);
        
    case 'helix'
        T = 10;
        t_vec = (0:dt:T)';
        
        om = 2*pi*0.125;
        x_d = sin(om*t_vec);
        y_d = cos(om*t_vec);
        z_d = -0.1*t_vec;
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = om*cos(om*t_vec);
        vy_d = -om*sin(om*t_vec);
        vz_d = -0.1*ones(length(t_vec),1);
        yd_d = zeros(length(t_vec),1);
        
        ax_d = -(om^2)*sin(om*t_vec);
        ay_d = -(om^2)*cos(om*t_vec);
        az_d = zeros(length(t_vec),1);
        
        jx_d = -(om^3)*cos(om*t_vec);
        jy_d =  (om^3)*sin(om*t_vec);
        jz_d = zeros(length(t_vec),1);
        
    case 'hover'
        T = 15;
        t_vec = (0:dt:T)';
        
        x_d = zeros(length(t_vec),1);
        y_d = zeros(length(t_vec),1);
        z_d = -0.5*ones(length(t_vec),1);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = zeros(length(t_vec),1);
        vy_d = zeros(length(t_vec),1);
        vz_d = zeros(length(t_vec),1);
        yd_d = zeros(length(t_vec),1);
        
        ax_d = zeros(length(t_vec),1);
        ay_d = zeros(length(t_vec),1);
        az_d = zeros(length(t_vec),1);
        
        jx_d = zeros(length(t_vec),1);
        jy_d = zeros(length(t_vec),1);
        jz_d = zeros(length(t_vec),1);
        
    case 'polyspline'
%         load('hardware_traj.mat');
        load(poly_file);
        p_order = 10;
        
        %Create splines for each segment
        slow_scale = 1.0;
        break_T = [0;cumsum(T_seg)];
        t_vec = (0:dt:break_T(end)*slow_scale)';
        x_d = zeros(length(t_vec),1); y_d = x_d; z_d = x_d;
        vx_d = x_d; vy_d = x_d; vz_d = x_d;
        ax_d = x_d; ay_d = x_d; az_d = x_d;
        jx_d = x_d; jy_d = x_d; jz_d = x_d;
        ns = 1;
        for i = 1:length(t_vec)
            t = t_vec(i)/slow_scale;
            if t > break_T(ns+1)
                ns = ns+1;
            end
            t_off = t - break_T(ns);
            pos_d = pos_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
            x_d(i) = pos_d(1); y_d(i) = pos_d(2); z_d(i) = pos_d(3);
            
            vel_d = (1/slow_scale)*vel_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
            vx_d(i) = vel_d(1); vy_d(i) = vel_d(2); vz_d(i) = vel_d(3);
            
            acc_d = (1/slow_scale^2)*acc_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
            ax_d(i) = acc_d(1); ay_d(i) = acc_d(2); az_d(i) = acc_d(3);
            
            jer_d = (1/slow_scale^3)*jer_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
            jx_d(i) = jer_d(1); jy_d(i) = jer_d(2); jz_d(i) = jer_d(3);
        end
        
        yaw_d = zeros(length(t_vec),1);
        yd_d = zeros(length(t_vec),1);
        
    otherwise
        disp('Unexpected traj type');        
end
%% Plot
figure()
plot3(x_d,y_d,z_d,'linewidth',2); grid on
xlabel('x'); ylabel('y'); zlabel('h');
set(gca,'ZDir','Reverse');
set(gca,'YDir','Reverse');

%% Mapping

%Body rates to Euler rates
R_eul = @(q) [cos(q(3))/cos(q(2)), -sin(q(3))/cos(q(2)), 0;
    sin(q(3)), cos(q(3)), 0;
    -cos(q(3))*tan(q(2)), sin(q(3))*tan(q(2)), 1];

%Get control & nominal state
u_nom = zeros(length(t_vec),4); %th_dot, om_x, om_y, om_z
x_nom = zeros(length(t_vec),9); %x,y,z, vx, vy,vz, r, p, y
th_nom = zeros(length(t_vec),1);

x_nom(:,1:6) = [x_d,y_d,z_d,vx_d,vy_d,vz_d];
x_nom(:,9) = yaw_d;

for i = 1:length(t_vec)
    th_vect = [ax_d(i);ay_d(i);az_d(i)-g];
    zb_ax = -th_vect/norm(th_vect);
    
    thrust = mq*norm(th_vect);
    thrust_d = -mq*(zb_ax'*[jx_d(i);jy_d(i);jz_d(i)]);
    th_nom(i) = thrust; %thrust
    u_nom(i,1) = thrust_d ; %thrust_dot
    
    roll = atan2(ay_d(i),g-az_d(i));
    pitch = asin(-ax_d(i)/(thrust/mq));
    
    x_nom(i,7) = roll;
    x_nom(i,8) = pitch;
    
    R = rot_matrix(roll,pitch,yaw_d(i));
    
    om_1 = [jx_d(i);jy_d(i);jz_d(i)]'*(R*[0;1;0])/(thrust/mq);
    om_2 = -[jx_d(i);jy_d(i);jz_d(i)]'*(R*[1;0;0])/(thrust/mq);
    
    roll_d = om_1*(cos(yaw_d(i))/cos(pitch)) - om_2*(sin(yaw_d(i))/cos(pitch));
    om_3 = yd_d(i) + roll_d*sin(pitch);
    
    om = [om_1;om_2;om_3];
    
    u_nom(i,2:4) = (R_eul([roll;pitch;yaw_d(i)])*om)';
    
end

accel = [ax_d,ay_d,az_d];

end