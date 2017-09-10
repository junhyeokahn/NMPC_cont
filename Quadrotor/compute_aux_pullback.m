function [uc_fb, E] = compute_aux_pullback(f_ctrl,B_ctrl,lambda, M_xi, xc_nom, yaw_nom, uc_nom, xc, yaw,uc_aux_prev)


%% Definitions

xic_nom = phi(xc_nom);
xic = phi(xc);
E = (xic-xic_nom)'*M_xi*(xic-xic_nom);

%% Compute control

d_xic = xic-xic_nom; %geodesic velocity in xi space

M_nom = M_fnc(xc_nom, M_xi);
M_act = M_fnc(xc, M_xi);

b = (2*d_xic'*(M_act*B_ctrl))';

a = 2*lambda*E - 2*d_xic'*(M_nom*(f_ctrl(xc_nom) + B_ctrl*uc_nom)) +...
                 2*d_xic'*(M_act*(f_ctrl(xc) + B_ctrl*(uc_nom+uc_aux_prev)));
                               
yd = 4*lambda*(yaw_nom-yaw); %simple proportional control or can get fancy

% a + b'aux <= 0
% a = (abs(a)>1e-4).*a;
% b = (abs(b)>1e-4).*b;

if (a <= 0) 
    u_fb_up = zeros(3,1); %[th_dot;rd;pd]
else
    u_fb_up = -(a/(b'*b))*b; %[th_dot;rd;pd]
end

u_metric_fb = uc_aux_prev + u_fb_up;

uc_fb = [u_metric_fb',yd];

end