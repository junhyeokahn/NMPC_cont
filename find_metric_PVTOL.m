function [solved,w_lower,w_upper,W_mat,dW_p_mat,dW_vy_mat,dW_vz_mat,dW_pd_mat,W_upper_mat] = ...
 find_metric_PVTOL(n,g,p_lim,pd_lim,vy_lim,vz_lim,...
    condn,lambda,ccm_eps,return_metric)
%% 


W_scale = diag([0.02;0.0268;0.0367;0.0089;0.023;0.02]);
% W_scale = 0.01*diag([20;10;5;1;20;15]);

% sin_x = @(x) 0.05059*(x/(pi/6));
% cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

% sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
% cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

yalmip('clear');

w_upper = sdpvar(1);
w_lower = sdpvar(1);

%states
Y = sdpvar(1);
Z = sdpvar(1);
p = sdpvar(1);
vy = sdpvar(1);
vz = sdpvar(1);
pd = sdpvar(1);

sin_p = sin_x(p);
cos_p = cos_x(p);

%dynamics f
f = [vy*cos_p - vz*sin_p;
     vy*sin_p + vz*cos_p;
     pd;
     pd*vz-g*sin_p;
     -pd*vy-g*cos_p;
     0];

%          Y Z  p                  vy    vz   pd
df_perp = [0,0,-vy*sin_p-vz*cos_p,cos_p,-sin_p,0;
           0,0, vy*cos_p-vz*sin_p,sin_p,cos_p,0;
           0,0,0,0,0,1;
           0,0,-g*cos_p,0,pd,vz];

B_perp = [eye(4);
          zeros(2,4)];

%Initialize decision variables
dec_List = [];

%% Parametrize W (2)

% w_states = [vy, p];

w_order = 2;
% W = sdpvar(6);
% for i = 1:n
%     for j = 1:n
%         if (j>=i)
%             [W(i,j),cp_c] = polynomial([p,vy],w_order);
%             dec_List = [dec_List;cp_c];
%         else
%             W(i,j) = W(j,i);
%         end
%     end
% end
W_perp = sdpvar(4);
for i = 1:4
    for j = 1:4
        if (j>=i)
            [W_perp(i,j),cp_c] = polynomial([p,vy],4);
            dec_List = [dec_List;cp_c];
        else
            W_perp(i,j) = W_perp(j,i);
        end
    end
end

W_pc = sdpvar(4,2);
for i = 1:4
    for j = 1:2
        [W_pc(i,j),cp_c] = polynomial([p,vy],w_order);
        dec_List = [dec_List;cp_c];
    end
end

W_c = sdpvar(2);
for i = 1:2
    for j = 1:2
        if (j>=i)
            [W_c(i,j),cp_c] = polynomial([p,vy],w_order);
            dec_List = [dec_List;cp_c];
        else
            W_c(i,j) = W_c(j,i);
        end
    end
end

W = [W_perp, W_pc;
     W_pc', W_c];
 
dW_f_perp = jacobian(W_perp(:),[Y,Z,p,vy,vz,pd])*f;
dW_f_perp = reshape(dW_f_perp,4,4);

W_upper = sdpvar(n);
dec_List = [dec_List; reshape(W_upper,n*n,1)];

%%

%Lagrange multipliers
% box_lim = [p+p_lim; p_lim-p;
%            vy+vy_lim; vy_lim-vy;
%            vz+vz_lim; vz_lim-vz;
%            pd+pd_lim; pd_lim-pd];
box_lim = [p_lim^2-p^2;
           vy_lim^2-vy^2;
           vz_lim^2-vz^2;
           pd_lim^2-pd^2];

delta_6 = sdpvar(6,1);
l_order = 4;
l_def_states = [p;vy; delta_6]';

[Ll_1, cl_1] = polynomial(l_def_states,l_order);
[Ll_2, cl_2] = polynomial(l_def_states,l_order);
% [Ll_3, cl_3] = polynomial(l_def_states,l_order);
% [Ll_4, cl_4] = polynomial(l_def_states,l_order);
Ll = [Ll_1,Ll_2];%,Ll_3,Ll_4];
dec_List = [dec_List;cl_1;cl_2];%cl_3;cl_4];

[Lu_1, cu_1] = polynomial(l_def_states,l_order);
[Lu_2, cu_2] = polynomial(l_def_states,l_order);
% [Lu_3, cu_3] = polynomial(l_def_states,l_order);
% [Lu_4, cu_4] = polynomial(l_def_states,l_order);
Lu = [Lu_1,Lu_2];%,Lu_3,Lu_4];
dec_List = [dec_List;cu_1;cu_2];%cu_3;cu_4];

delta_4 = sdpvar(4,1);
l_ccm_states = [p;vy;pd;vz;delta_4]';
lc_order = 4;

[Lc_1, cc_1] = polynomial(l_ccm_states,lc_order+2);
[Lc_2, cc_2] = polynomial(l_ccm_states,lc_order);
[Lc_3, cc_3] = polynomial(l_ccm_states,lc_order);
[Lc_4, cc_4] = polynomial(l_ccm_states,lc_order);
% [Lc_5, cc_5] = polynomial(l_ccm_states,lc_order);
% [Lc_6, cc_6] = polynomial(l_ccm_states,lc_order);
% [Lc_7, cc_7] = polynomial(l_ccm_states,lc_order);
% [Lc_8, cc_8] = polynomial(l_ccm_states,lc_order);
Lc = [Lc_1,Lc_2,Lc_3,Lc_4];%,Lc_5,Lc_6,Lc_7,Lc_8];
dec_List = [dec_List;cc_1;cc_2;cc_3;cc_4];%cc_5;cc_6;cc_7;cc_8];

%W uniform bounds
W_bounds = [w_lower>=1.0, W_upper <= w_upper*eye(n)];

%Condition bound
W_cond = [w_upper-condn*w_lower <= 0];

%W pos def
p_low = (delta_6'*W*delta_6 - w_lower*(delta_6'*delta_6)) - Ll*box_lim(1:2);
p_high = delta_6'*(W_upper - W)*delta_6 - Lu*box_lim(1:2);

%CCM condition
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);
p_CCM = (delta_4'*R_CCM*delta_4 - ccm_eps*(delta_4'*delta_4)) - Lc*box_lim;

%Assemble
constr_List = [W_bounds; W_cond; sos(p_low); sos(p_high); sos(p_CCM);
    sos(Ll_1); sos(Ll_2);% sos(Ll_3); sos(Ll_4);
    sos(Lu_1); sos(Lu_2);% sos(Lu_3); sos(Lu_4);
    sos(Lc_1); sos(Lc_2); sos(Lc_3); sos(Lc_4)];

options = sdpsettings('solver','mosek','verbose',return_metric);
SOS_soln = solvesos(constr_List,trace(W_scale*W_upper) + (1e-3)*norm(dec_List,1),options, dec_List);

%log(det(W_scale'*W_upper*W_scale))
solved = SOS_soln.problem;

%% Parse
w_lower = double(w_lower);
w_upper = double(w_upper);

W_mat = zeros(n);
dW_p_mat = zeros(n);
dW_vy_mat = zeros(n);

dW_vz_mat = zeros(n);
dW_pd_mat = zeros(n);

W_upper_mat = zeros(n);

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');
        dec_sol = clean(double(dec_List),1e-3);
        W_sol = replace(W,dec_List,dec_sol);
        
        dW_p = jacobian(W(:),p);
        dW_p_sol = replace(reshape(dW_p,n,n),dec_List,dec_sol);
        
        dW_vy = jacobian(W(:),vy);
        dW_vy_sol = replace(reshape(dW_vy,n,n),dec_List,dec_sol);       
        
%         dW_vz = jacobian(W(:),vz);
%         dW_vz_sol = replace(reshape(dW_vz,n,n),dec_List,dec_sol);       
%         
%         dW_pd = jacobian(W(:),pd);
%         dW_pd_sol = replace(reshape(dW_pd,n,n),dec_List,dec_sol);       
        
        W_upper_mat = replace(W_upper,dec_List,dec_sol);
%         W_upper_mat = double(W_upper);
        sdisplay(W_sol)
        sdisplay(W_upper_mat);
        pause;
                
        %% Create functions

        W_cell = cell(n,n);
        dW_vy_cell = cell(n,n);
        dW_p_cell = cell(n,n);
        
%         dW_vz_cell = cell(n,n);
%         dW_pd_cell = cell(n,n);
        
        for i = 1:n
            for j = 1:n
                W_ij = sdisplay(W_sol(i,j));
                dW_vy_ij = sdisplay(dW_vy_sol(i,j));
                dW_p_ij = sdisplay(dW_p_sol(i,j));
%                 dW_vz_ij = sdisplay(dW_vz_sol(i,j));
%                 dW_pd_ij = sdisplay(dW_pd_sol(i,j));
                
                W_ij = strcat('@(p,vy)',W_ij{1});
                dW_vy_ij = strcat('@(p,vy)',dW_vy_ij{1});
                dW_p_ij = strcat('@(p,vy)',dW_p_ij{1});
%                 dW_vz_ij = strcat('@(p,vy,vz,pd)',dW_vz_ij{1});
%                 dW_pd_ij = strcat('@(p,vy,vz,pd)',dW_pd_ij{1});
                
                W_cell{i,j} = str2func(W_ij);
                dW_vy_cell{i,j} = str2func(dW_vy_ij);
                dW_p_cell{i,j} = str2func(dW_p_ij);
%                 dW_vz_cell{i,j} = str2func(dW_vz_ij);
%                 dW_pd_cell{i,j} = str2func(dW_pd_ij);
            end
        end
        W_mat = @(x) [W_cell{1,1}(x(3),x(4)), W_cell{1,2}(x(3),x(4)), W_cell{1,3}(x(3),x(4)), W_cell{1,4}(x(3),x(4)), W_cell{1,5}(x(3),x(4)), W_cell{1,6}(x(3),x(4));
                      W_cell{2,1}(x(3),x(4)), W_cell{2,2}(x(3),x(4)), W_cell{2,3}(x(3),x(4)), W_cell{2,4}(x(3),x(4)), W_cell{2,5}(x(3),x(4)), W_cell{2,6}(x(3),x(4));
                      W_cell{3,1}(x(3),x(4)), W_cell{3,2}(x(3),x(4)), W_cell{3,3}(x(3),x(4)), W_cell{3,4}(x(3),x(4)), W_cell{3,5}(x(3),x(4)), W_cell{3,6}(x(3),x(4));
                      W_cell{4,1}(x(3),x(4)), W_cell{4,2}(x(3),x(4)), W_cell{4,3}(x(3),x(4)), W_cell{4,4}(x(3),x(4)), W_cell{4,5}(x(3),x(4)), W_cell{4,6}(x(3),x(4));
                      W_cell{5,1}(x(3),x(4)), W_cell{5,2}(x(3),x(4)), W_cell{5,3}(x(3),x(4)), W_cell{5,4}(x(3),x(4)), W_cell{5,5}(x(3),x(4)), W_cell{5,6}(x(3),x(4));
                      W_cell{6,1}(x(3),x(4)), W_cell{6,2}(x(3),x(4)), W_cell{6,3}(x(3),x(4)), W_cell{6,4}(x(3),x(4)), W_cell{6,5}(x(3),x(4)), W_cell{6,6}(x(3),x(4))];
        
        dW_vy_mat = @(x) [dW_vy_cell{1,1}(x(3),x(4)), dW_vy_cell{1,2}(x(3),x(4)), dW_vy_cell{1,3}(x(3),x(4)), dW_vy_cell{1,4}(x(3),x(4)), dW_vy_cell{1,5}(x(3),x(4)), dW_vy_cell{1,6}(x(3),x(4));
                          dW_vy_cell{2,1}(x(3),x(4)), dW_vy_cell{2,2}(x(3),x(4)), dW_vy_cell{2,3}(x(3),x(4)), dW_vy_cell{2,4}(x(3),x(4)), dW_vy_cell{2,5}(x(3),x(4)), dW_vy_cell{2,6}(x(3),x(4));
                          dW_vy_cell{3,1}(x(3),x(4)), dW_vy_cell{3,2}(x(3),x(4)), dW_vy_cell{3,3}(x(3),x(4)), dW_vy_cell{3,4}(x(3),x(4)), dW_vy_cell{3,5}(x(3),x(4)), dW_vy_cell{3,6}(x(3),x(4));
                          dW_vy_cell{4,1}(x(3),x(4)), dW_vy_cell{4,2}(x(3),x(4)), dW_vy_cell{4,3}(x(3),x(4)), dW_vy_cell{4,4}(x(3),x(4)), dW_vy_cell{4,5}(x(3),x(4)), dW_vy_cell{4,6}(x(3),x(4));
                          dW_vy_cell{5,1}(x(3),x(4)), dW_vy_cell{5,2}(x(3),x(4)), dW_vy_cell{5,3}(x(3),x(4)), dW_vy_cell{5,4}(x(3),x(4)), dW_vy_cell{5,5}(x(3),x(4)), dW_vy_cell{5,6}(x(3),x(4));
                          dW_vy_cell{6,1}(x(3),x(4)), dW_vy_cell{6,2}(x(3),x(4)), dW_vy_cell{6,3}(x(3),x(4)), dW_vy_cell{6,4}(x(3),x(4)), dW_vy_cell{6,5}(x(3),x(4)), dW_vy_cell{6,6}(x(3),x(4))];
        
        dW_p_mat = @(x)  [dW_p_cell{1,1}(x(3),x(4)), dW_p_cell{1,2}(x(3),x(4)), dW_p_cell{1,3}(x(3),x(4)), dW_p_cell{1,4}(x(3),x(4)), dW_p_cell{1,5}(x(3),x(4)), dW_p_cell{1,6}(x(3),x(4));
                          dW_p_cell{2,1}(x(3),x(4)), dW_p_cell{2,2}(x(3),x(4)), dW_p_cell{2,3}(x(3),x(4)), dW_p_cell{2,4}(x(3),x(4)), dW_p_cell{2,5}(x(3),x(4)), dW_p_cell{2,6}(x(3),x(4));
                          dW_p_cell{3,1}(x(3),x(4)), dW_p_cell{3,2}(x(3),x(4)), dW_p_cell{3,3}(x(3),x(4)), dW_p_cell{3,4}(x(3),x(4)), dW_p_cell{3,5}(x(3),x(4)), dW_p_cell{3,6}(x(3),x(4));
                          dW_p_cell{4,1}(x(3),x(4)), dW_p_cell{4,2}(x(3),x(4)), dW_p_cell{4,3}(x(3),x(4)), dW_p_cell{4,4}(x(3),x(4)), dW_p_cell{4,5}(x(3),x(4)), dW_p_cell{4,6}(x(3),x(4));
                          dW_p_cell{5,1}(x(3),x(4)), dW_p_cell{5,2}(x(3),x(4)), dW_p_cell{5,3}(x(3),x(4)), dW_p_cell{5,4}(x(3),x(4)), dW_p_cell{5,5}(x(3),x(4)), dW_p_cell{5,6}(x(3),x(4));
                          dW_p_cell{6,1}(x(3),x(4)), dW_p_cell{6,2}(x(3),x(4)), dW_p_cell{6,3}(x(3),x(4)), dW_p_cell{6,4}(x(3),x(4)), dW_p_cell{6,5}(x(3),x(4)), dW_p_cell{6,6}(x(3),x(4))];
        
%         dW_vz_mat = @(x) [dW_vz_cell{1,1}(x(3),x(4)), dW_vz_cell{1,2}(x(3),x(4)), dW_vz_cell{1,3}(x(3),x(4)), dW_vz_cell{1,4}(x(3),x(4)), dW_vz_cell{1,5}(x(3),x(4)), dW_vz_cell{1,6}(x(3),x(4));
%                           dW_vz_cell{2,1}(x(3),x(4)), dW_vz_cell{2,2}(x(3),x(4)), dW_vz_cell{2,3}(x(3),x(4)), dW_vz_cell{2,4}(x(3),x(4)), dW_vz_cell{2,5}(x(3),x(4)), dW_vz_cell{2,6}(x(3),x(4));
%                           dW_vz_cell{3,1}(x(3),x(4)), dW_vz_cell{3,2}(x(3),x(4)), dW_vz_cell{3,3}(x(3),x(4)), dW_vz_cell{3,4}(x(3),x(4)), dW_vz_cell{3,5}(x(3),x(4)), dW_vz_cell{3,6}(x(3),x(4));
%                           dW_vz_cell{4,1}(x(3),x(4)), dW_vz_cell{4,2}(x(3),x(4)), dW_vz_cell{4,3}(x(3),x(4)), dW_vz_cell{4,4}(x(3),x(4)), dW_vz_cell{4,5}(x(3),x(4)), dW_vz_cell{4,6}(x(3),x(4));
%                           dW_vz_cell{5,1}(x(3),x(4)), dW_vz_cell{5,2}(x(3),x(4)), dW_vz_cell{5,3}(x(3),x(4)), dW_vz_cell{5,4}(x(3),x(4)), dW_vz_cell{5,5}(x(3),x(4)), dW_vz_cell{5,6}(x(3),x(4));
%                           dW_vz_cell{6,1}(x(3),x(4)), dW_vz_cell{6,2}(x(3),x(4)), dW_vz_cell{6,3}(x(3),x(4)), dW_vz_cell{6,4}(x(3),x(4)), dW_vz_cell{6,5}(x(3),x(4)), dW_vz_cell{6,6}(x(3),x(4))];
%         
%         dW_pd_mat = @(x) [dW_pd_cell{1,1}(x(3),x(4)), dW_pd_cell{1,2}(x(3),x(4)), dW_pd_cell{1,3}(x(3),x(4)), dW_pd_cell{1,4}(x(3),x(4)), dW_pd_cell{1,5}(x(3),x(4)), dW_pd_cell{1,6}(x(3),x(4));
%                           dW_pd_cell{2,1}(x(3),x(4)), dW_pd_cell{2,2}(x(3),x(4)), dW_pd_cell{2,3}(x(3),x(4)), dW_pd_cell{2,4}(x(3),x(4)), dW_pd_cell{2,5}(x(3),x(4)), dW_pd_cell{2,6}(x(3),x(4));
%                           dW_pd_cell{3,1}(x(3),x(4)), dW_pd_cell{3,2}(x(3),x(4)), dW_pd_cell{3,3}(x(3),x(4)), dW_pd_cell{3,4}(x(3),x(4)), dW_pd_cell{3,5}(x(3),x(4)), dW_pd_cell{3,6}(x(3),x(4));
%                           dW_pd_cell{4,1}(x(3),x(4)), dW_pd_cell{4,2}(x(3),x(4)), dW_pd_cell{4,3}(x(3),x(4)), dW_pd_cell{4,4}(x(3),x(4)), dW_pd_cell{4,5}(x(3),x(4)), dW_pd_cell{4,6}(x(3),x(4));
%                           dW_pd_cell{5,1}(x(3),x(4)), dW_pd_cell{5,2}(x(3),x(4)), dW_pd_cell{5,3}(x(3),x(4)), dW_pd_cell{5,4}(x(3),x(4)), dW_pd_cell{5,5}(x(3),x(4)), dW_pd_cell{5,6}(x(3),x(4));
%                           dW_pd_cell{6,1}(x(3),x(4)), dW_pd_cell{6,2}(x(3),x(4)), dW_pd_cell{6,3}(x(3),x(4)), dW_pd_cell{6,4}(x(3),x(4)), dW_pd_cell{6,5}(x(3),x(4)), dW_pd_cell{6,6}(x(3),x(4))];
        
    end
end
end