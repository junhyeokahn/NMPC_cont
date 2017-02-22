function [solved,w_lower,w_upper] = ...
    find_metric_QUAD_SPOT_89_c(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,condn,lambda,ccm_eps,return_metric)
%%


% W_scale = diag([0.02;0.0268;0.0367;0.0089;0.02;0.02]);
% W_scale = diag([0.001;0.02;0.0367;0.005;0.001;0.002]);
% W_scale = zeros(n);

sinr_lim = sin(r_lim);
cosr_lim_low = cos(r_lim);
sinp_lim = sin(p_lim);
cosp_lim_low = cos(p_lim);

%states
n = 9;
xt = msspoly('xt',7); %pos, vel, thrust
sinr = msspoly('sinr',1);
cosr = msspoly('cosr',1);
sinp = msspoly('sinp',1);
cosp = msspoly('cosp',1);

%pos_def indeterminates
dnin = msspoly('dnin',9);
dsix = msspoly('dsix',6);

%dynamics f (xt)
b_T = [sinp; -cosp*sinr; cosp*cosr];

f = [[xt(4);xt(5);xt(6)];
     [0;0;g] - b_T*xt(7);
      0]; 

%gradients       
db_T_q = [0, cosp;
         -cosr*cosp, sinr*sinp;
         -sinr*cosp,-cosr*sinp];

%          x y z vx     vy    vz t r p       
df_perp      = [zeros(3), eye(3),zeros(3,3);
                zeros(3,6), -b_T, -db_T_q(:,1)*xt(7), -db_T_q(:,2)*xt(7)];
            
B_perp = [eye(6);
          zeros(3,6)];
      
%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(xt);
prog = prog.withIndeterminate(sinr);
prog = prog.withIndeterminate(cosr);
prog = prog.withIndeterminate(sinp);
prog = prog.withIndeterminate(cosp);
prog = prog.withIndeterminate(dnin);
prog = prog.withIndeterminate(dsix);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);

%% Parametrize W (2)

w_states = [xt(7);sinr;cosr;sinp;cosp];

w_order = 2;
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

[prog, W_upper] = prog.newSym(n);

%% Killing Field conditions

% kill_1 = diff(W_perp(:),xt(7));
% kill_2 = diff(W_perp(:),sinr)*cosr + diff(W_perp(:),cosr)*(-sinr);
% kill_3 = diff(W_perp(:),sinp)*cosp + diff(W_perp(:),cosp)*(-sinp);

% prog = prog.withPolyEqs([kill_1(:);
%                          kill_2(:);
%                          kill_3(:)]);

%% Definiteness conditions

%Lagrange multipliers
box_lim = [%vx_lim^2-xt(4)^2;
           %vy_lim^2-xt(5)^2;
           %vz_lim^2-xt(6)^2;
           sinr_lim^2-sinr^2;
           sinp_lim^2-sinp^2;
           cosr-cosr_lim_low;
           cosp-cosp_lim_low;
           xt(7) - th_lim_low;
           th_lim_high - xt(7)];
       
sphere_def = [cosr^2 + sinr^2 - 1;
              cosp^2 + sinp^2 - 1];

%For definiteness          
l_order = 2;
l_def_states = [w_states; dnin]';
n_def_L = length(box_lim);

[prog, Ll] = prog.newSDSOSPoly(monomials(l_def_states,0:l_order),n_def_L);
[prog, Lu] = prog.newSDSOSPoly(monomials(l_def_states,0:l_order),n_def_L);
[prog, Ll_sph] = prog.newFreePoly(monomials(l_def_states,0:l_order),2);
[prog, Lu_sph] = prog.newFreePoly(monomials(l_def_states,0:l_order),2);

%For CCM
lc_order = 2;
l_ccm_states = [w_states;dsix]';

% [prog, Lc_v]  = prog.newSDSOSPoly(monomials(l_ccm_states,0:lc_order),3);
[prog, Lc_rp] = prog.newSDSOSPoly(monomials(l_ccm_states,0:lc_order+2),4);
[prog, Lc_th] = prog.newSDSOSPoly(monomials(l_ccm_states,0:lc_order),2);
Lc = [Lc_rp; Lc_th];
[prog, Lc_sph] = prog.newFreePoly(monomials(l_ccm_states,0:lc_order+2),2);

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
% prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSDSOS((dnin'*W*dnin - w_lower*(dnin'*dnin)) - Ll'*box_lim(1:n_def_L) + Ll_sph'*sphere_def);
prog = prog.withSDSOS(dnin'*(W_upper - W)*dnin - Lu'*box_lim(1:n_def_L)+ Lu_sph'*sphere_def);

%CCM condition
R_CCM = -(-dW_perp_f + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);
prog = prog.withSOS((dsix'*R_CCM*dsix - ccm_eps*(dsix'*dsix)) - Lc'*box_lim + Lc_sph'*sphere_def);

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.coneVar; prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

% SOS_soln = prog.minimize(trace(W_scale*W_upper) + (1e-3)*sum(a), @spot_mosek, options);
SOS_soln = prog.minimize(w_upper-w_lower + (1e-4)*sum(a), @spot_mosek, options);

solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

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
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-3);
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
        
        W_upper_mat = clean(double(SOS_soln.eval(W_upper)),1e-3);
        
        pause;
        
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
        save('metric_QUAD_vectorized.mat','W_eval','w_poly_fnc','dw_poly_vx_fnc','dw_poly_vy_fnc','dw_poly_vz_fnc','dw_poly_r_fnc','dw_poly_p_fnc','dw_poly_th_fnc')
%         eval(W_exec); eval(dW_x1_exec); eval(dW_x2_exec);
    end
end
end