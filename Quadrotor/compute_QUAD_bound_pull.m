function [m_lower, M_lower_pull] = compute_QUAD_bound_pull(M_xi,r_lim,p_lim,th_lim_low,th_lim_high)

% sin_x = @(x) 0.5059*(x/(pi/6));
% cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

% sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

dnin = msspoly('dnin',9);
r = msspoly('r',1);
p = msspoly('p',1);
th = msspoly('th',1);

sin_r = sin_x(r);
cos_r = cos_x(r);
sin_p = sin_x(p);
cos_p = cos_x(p);

b_T = [sin_p; -cos_p*sin_r; cos_p*cos_r];
db_T_q = [0, cos_p;
         -cos_r*cos_p, sin_r*sin_p;
         -sin_r*cos_p,-cos_r*sin_p];
Phi =  blkdiag(eye(3),eye(3),-[b_T, db_T_q*th]);

M = Phi'*M_xi*Phi;

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(r);
prog = prog.withIndeterminate(p);
prog = prog.withIndeterminate(th);
prog = prog.withIndeterminate(dnin);

[prog, m_lower] = prog.newPos(1);

%% Variable

[prog, M_lower_pull] = prog.newSym(9);

%% Lagrange multipliers

box_lim = [r_lim^2-r^2;
           p_lim^2-p^2;
           th - th_lim_low;
           th_lim_high - th];
       
l_order = 4;
l_def_states = [r;p;th];
n_def_L = length(box_lim);

[m_def_mon, m_def_mat] = monomials([l_def_states;dnin],0:l_order);
m_def_keep = find(sum(m_def_mat(:,length(l_def_states)+1:end),2)==2); %only keep quadratics in dnin
m_def_mon = m_def_mon(m_def_keep);

[prog, Ll] = prog.newSOSPoly(m_def_mon,n_def_L);

%Bounds
prog = prog.withPSD(M_lower_pull - m_lower*eye(9));
prog = prog.withSOS(dnin'*(M - M_lower_pull)*dnin - (Ll'*box_lim));

options = spot_sdp_default_options();
options.verbose = 1;

SOS_soln = prog.minimize(-m_lower - 0.1*trace(M_lower_pull), @spot_mosek, options);
% 0.1*trace(M_lower_pull)
solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

if (solved == 0)
    m_lower = double(SOS_soln.eval(m_lower));
    M_lower_pull = clean(double(SOS_soln.eval(M_lower_pull)),1e-6);    
end
    
end
