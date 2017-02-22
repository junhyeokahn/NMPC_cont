function [solved,w_lower,w_upper] = ...
    find_metric_QUAD_SPOT_89(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,condn,lambda,ccm_eps,return_metric)
%%


W_scale = (1e-1)*diag([0.002*ones(6,1);2;2;0.005]);
% W_scale = zeros(9);
norm_scale = 1e-3;

load 'quad_metric_pull.mat'
% L = chol(W_upper_pull);

sin_x = @(x) 0.5059*(x/(pi/6));
cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

% sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

% sin_yaw = @(x) 0.5692*(x/pi) - 0.6669*(4*(x/pi)^3 - 3*(x/pi)) + ...
%     0.1043*(16*(x/pi)^5 - 20*(x/pi)^3 + 5*(x/pi));

% cos_yaw = @(x)   -0.3042 - 0.9709*(2*(x/pi)^2 -1) + 0.3028*(8*(x/pi)^4 - 8*(x/pi)^2 + 1) + ...
%     -0.029*(32*(x/pi)^6 - 48*(x/pi)^4 + 18*(x/pi)^2 - 1);

% sin_yaw = sin_x;
% cos_yaw = cos_x;

%states
n = 9;
x = msspoly('x',9);

%pos_def indeterminates
dnin = msspoly('dnin',9);
dsix = msspoly('dsix',6);

sin_r = sin_x(x(7));
cos_r = cos_x(x(7));
sin_p = sin_x(x(8));
cos_p = cos_x(x(8));

%dynamics f
b_T = [sin_p; -cos_p*sin_r; cos_p*cos_r];

f = [[x(4);x(5);x(6)];
     [0;0;g] - b_T*x(9);
     zeros(3,1)]; 

%gradients       
db_T_q = [0, cos_p;
         -cos_r*cos_p, sin_r*sin_p;
         -sin_r*cos_p,-cos_r*sin_p];

%          x y z vx     vy    vz  r p  t     
df_perp      = [zeros(3), eye(3),zeros(3,3);
                zeros(3,6), -db_T_q(:,1)*x(9), -db_T_q(:,2)*x(9), -b_T];
            
B_perp = [eye(6);
          zeros(3,6)];
      
%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dnin);
prog = prog.withIndeterminate(dsix);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);
[prog, W_upper] = prog.newSym(n);

%% Parametrize W (2)

w_states = [x(7);x(8);x(9)];

w_order = 4;
w_poly = monomials(w_states,0:w_order);
W_list = cell(length(w_poly),1);
W_pc_list = cell(length(w_poly),1);
W_c_list = cell(length(w_poly),1);

[prog, W_perp] = prog.newSym(6);
[prog, W_pc_list{1}] = prog.newFree(6,3);
[prog, W_c_list{1}] = prog.newSym(3);

W_list{1} = [W_perp, W_pc_list{1};
            W_pc_list{1}', W_c_list{1}];
W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    [prog, W_c_list{i}] = prog.newSym(3);
    [prog, W_pc_list{i}] = prog.newFree(6,3);
    W_list{i} = [zeros(6), W_pc_list{i};
                W_pc_list{i}', W_c_list{i}];
    W = W + W_list{i}*w_poly(i);
end

% W_perp = W(1:6,1:6);
% dW_perp_f = diff(W_perp(:),x)*f;
% dW_perp_f = reshape(dW_perp_f,6,6);
dW_perp_f = zeros(6);

%% Killing Field conditions

% kill_1 = diff(W_perp(:),x(7));
% kill_2 = diff(W_perp(:),x(8));
% kill_3 = diff(W_perp(:),x(9));
%  
% prog = prog.withPolyEqs([kill_1(:);
%                          kill_2(:);
%                          kill_3(:)]);

%% Definiteness conditions

%Lagrange multipliers
box_lim = [%vx_lim^2-x(4)^2;
           %vy_lim^2-x(5)^2;
           %vz_lim^2-x(6)^2;
           r_lim^2-x(7)^2;
           p_lim^2-x(8)^2;
           x(9) - th_lim_low;
           th_lim_high - x(9)];

l_order = 4;
l_def_states = [w_states; dnin]';
n_def_L = 4;

[prog, Ll] = prog.newSDSOSPoly(monomials(l_def_states,0:l_order),n_def_L);
[prog, Lu] = prog.newSDSOSPoly(monomials(l_def_states,0:l_order),n_def_L);

l_ccm_states = [w_states;dsix]';
lc_order = 4;

% [prog, Lc_v]  = prog.newSDSOSPoly(monomials(l_ccm_states,0:2),3);
[prog, Lc_rp] = prog.newSDSOSPoly(monomials(l_ccm_states,0:lc_order+2),2);
[prog, Lc_th] = prog.newSDSOSPoly(monomials(l_ccm_states,0:lc_order),2);
Lc = [Lc_rp; Lc_th];

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSDSOS((dnin'*W*dnin - w_lower*(dnin'*dnin)) - Ll'*box_lim(1:n_def_L));
prog = prog.withSDSOS(dnin'*(W_upper - W)*dnin - Lu'*box_lim(1:n_def_L));

%CCM condition
R_CCM = -(-dW_perp_f + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);
prog = prog.withSDSOS((dsix'*R_CCM*dsix - ccm_eps*(dsix'*dsix)) - Lc'*box_lim);

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.freeVar(0.5*n*(n+1)+1:end)];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

try
    SOS_soln = prog.minimize(trace(W_scale*W_upper) + norm_scale*sum(a), @spot_mosek, options);
catch
    %failed
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end

solved = ~(strcmp(SOS_soln.info.solverInfo.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE'));% && ...
%            strcmp(SOS_soln.info.solverInfo.itr.solsta, 'OPTIMAL'));

%% Parse

w_lower = double(SOS_soln.eval(w_lower));
w_upper = double(SOS_soln.eval(w_upper));

W_upper_mat = zeros(n);

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));
        NNZ_list = zeros(length(w_poly),1);
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-4);
            if sum(sum(abs(W_sol(:,:,i)))) > 0
                NNZ_list(i) = 1;
            end
        end
        w_poly = w_poly(find(NNZ_list));
        W_sol = W_sol(:,:,find(NNZ_list));
        
        dw_poly_vx = diff(w_poly,x(4));
        dw_poly_vy = diff(w_poly,x(5));
        dw_poly_vz = diff(w_poly,x(6));
        dw_poly_r = diff(w_poly,x(7));
        dw_poly_p = diff(w_poly,x(8));
        dw_poly_th = diff(w_poly,x(9));
        
        W_upper_mat = clean(double(SOS_soln.eval(W_upper)),1e-4);
        
%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_vx_fnc = mss2fnc(dw_poly_vx,x,randn(length(x),2));
        dw_poly_vy_fnc = mss2fnc(dw_poly_vy,x,randn(length(x),2));
        dw_poly_vz_fnc = mss2fnc(dw_poly_vz,x,randn(length(x),2));
        dw_poly_r_fnc = mss2fnc(dw_poly_r,x,randn(length(x),2));
        dw_poly_p_fnc = mss2fnc(dw_poly_p,x,randn(length(x),2));
        dw_poly_th_fnc = mss2fnc(dw_poly_th,x,randn(length(x),2));
        
        %% Put together
        W_exec = 'W_eval = @(ml)';
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute
        eval(W_exec);
        save('metric_QUAD_vectorized.mat','W_eval','w_poly_fnc','dw_poly_vx_fnc','dw_poly_vy_fnc','dw_poly_vz_fnc','dw_poly_r_fnc','dw_poly_p_fnc','dw_poly_th_fnc','W_upper_mat');
%         eval(W_exec); eval(dW_x1_exec); eval(dW_x2_exec);
    end
end
end