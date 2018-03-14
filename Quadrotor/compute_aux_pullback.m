function [uc_fb, E] = compute_aux_pullback(f_ctrl,B_ctrl,lambda, M_xi, xc_nom, yaw_nom, u_nom, xc, yaw)


%% Definitions
%ignore roll pitch errors below 0.1 deg
% rot_eps = 0.1*(pi/180);
% for i = 8:9
%     if (abs(xc(i)-xc_nom(i))<rot_eps)
%         xc(i) = xc_nom(i);
%     end
% end
xic_nom = phi(xc_nom);
xic = phi(xc);
uc_nom = u_nom(1:3);
E = (xic-xic_nom)'*M_xi*(xic-xic_nom);

%% Compute control

d_xic = xic-xic_nom; %geodesic velocity in xi space

M_nom = M_fnc(xc_nom, M_xi);
M_act = M_fnc(xc, M_xi);

b = (2*d_xic'*(M_act*B_ctrl))';

a = 2*lambda*E - 2*d_xic'*(M_nom*(f_ctrl(xc_nom) + B_ctrl*uc_nom)) +...
                 2*d_xic'*(M_act*(f_ctrl(xc) + B_ctrl*uc_nom));

% a + b'aux <= 0

if (a <= 0) 
    uc_fb = zeros(3,1); %[th_dot;rd;pd]
else
    uc_fb = -(a/(b'*b))*b; %[th_dot;rd;pd]
end

uc_fb = [uc_fb;-2*(yaw-yaw_nom)];

% uc_fb = uc_aux_prev + u_fb_up;

end