function [solved, w_lower, w_upper, W_mat, dW_x1_mat, dW_x2_mat,W_upper_mat] = find_metric_FLR(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
    condn, lambda, ccm_eps,return_metric)


W_scale = diag([0.01; 0.0005; 0.01; 0.0005]);
% W_scale = diag(zeros(4,1));

sin_x = @(x) 0.5692*(x/pi) - 0.6669*(4*(x/pi)^3 - 3*(x/pi)) + ...
    0.1043*(16*(x/pi)^5 - 20*(x/pi)^3 + 5*(x/pi));

cos_x = @(x)   -0.3042 - 0.9709*(2*(x/pi)^2 -1) + 0.3028*(8*(x/pi)^4 - 8*(x/pi)^2 + 1) + ...
    -0.029*(32*(x/pi)^6 - 48*(x/pi)^4 + 18*(x/pi)^2 - 1);

yalmip('clear');

w_upper = sdpvar(1);
w_lower = sdpvar(1);

%states
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);

sin_x1 = sin_x(x1);
cos_x1 = cos_x(x1);

%dynamics f
f = [x2;
    (m*g*l/I)*sin_x1 - (sigma/I)*(x1-x3);
    x4;
    (sigma/J)*(x1-x3) - (b/J)*x4];

df_perp = [0, 1, 0, 0;
    (m*g*l/I)*cos_x1-(sigma/I), 0, (sigma/I), 0;
    zeros(1,3), 1];

B_perp = [eye(3);
    zeros(1,3)];

%Initialize decision variables
dec_List = [];

%% Parametrize W (2)

% w_states = [x1,x2];
w_states = [x1];
w_order = 2;

w_mono = monolist(w_states,w_order);
W_list = cell(length(w_mono),1);

W_list{1} = sdpvar(n);
dec_List = [dec_List;reshape(W_list{1},n*n,1)];
W = W_list{1}*w_mono(1);

for i = 2:length(w_mono)
    W_list{i} = sdpvar(n);
    dec_List = [dec_List;reshape(W_list{i},n*n,1)];
    W = W + W_list{i}*w_mono(i);
end

% W = sdpvar(n);
% for i = 1:n
%     for j = 1:n
%         if (j>=i)
%             [W(i,j),cp_c] = polynomial(w_states,w_order);
%             dec_List = [dec_List;cp_c];
%         else
%             W(i,j) = W(j,i);
%         end
%     end
% end

W_upper = sdpvar(n);
dec_List = [dec_List; reshape(W_upper,n*n,1)];

dW_f = jacobian(W(:),[x1,x2,x3,x4])*f;
dW_f = reshape(dW_f,n,n);

%%

%Lagrange multipliers
box_lim = [x1+x1_lim; x1_lim-x1;
    x2+x2_lim; x2_lim-x2];
%            x3+x3_lim; x3_lim-x3];

delta_4 = sdpvar(4,1);
l_order = 4;
l_def_states = [w_states'; delta_4]';

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

delta_3 = sdpvar(3,1);
l_ccm_states = [x1;x2;delta_3]';
lc_order = 6;

[Lc_1, cc_1] = polynomial(l_ccm_states,lc_order);
[Lc_2, cc_2] = polynomial(l_ccm_states,lc_order);
[Lc_3, cc_3] = polynomial(l_ccm_states,lc_order);
[Lc_4, cc_4] = polynomial(l_ccm_states,lc_order);
% [Lc_5, cc_5] = polynomial(l_ccm_states,lc_order);
% [Lc_6, cc_6] = polynomial(l_ccm_states,lc_order);
Lc = [Lc_1,Lc_2,Lc_3,Lc_4];%,Lc_5,Lc_6];
dec_List = [dec_List;cc_1;cc_2;cc_3;cc_4];%cc_5;cc_6];

%W uniform bounds
W_bounds = [w_lower>=1, W_upper <= w_upper*eye(n)];

%Condition bound
W_cond = [w_upper-condn*w_lower <= 0];

%W pos def
p_low = (delta_4'*W*delta_4 - w_lower*(delta_4'*delta_4)) - Ll*box_lim(1:2);
p_high = delta_4'*(W_upper - W)*delta_4 - Lu*box_lim(1:2);

%CCM condition
R_CCM = -(-dW_f(1:3,1:3) + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W(1:3,1:3));
p_CCM = (delta_3'*R_CCM*delta_3 - ccm_eps*(delta_3'*delta_3)) - Lc*box_lim;

%Assemble
constr_List = [W_bounds; W_cond; sos(p_low); sos(p_high); sos(p_CCM);
    sos(Ll_1); sos(Ll_2);% sos(Ll_3); sos(Ll_4);
    sos(Lu_1); sos(Lu_2);% sos(Lu_3); sos(Lu_4);
    sos(Lc_1); sos(Lc_2); sos(Lc_3); sos(Lc_4)];% sos(Lc_5); sos(Lc_6)];

options = sdpsettings('solver','mosek','verbose',return_metric);
SOS_soln = solvesos(constr_List,trace(W_scale*W_upper) + (1e-3)*norm(dec_List,1),options, dec_List);

solved = SOS_soln.problem;
w_lower = double(w_lower);
w_upper = double(w_upper);

W_mat = zeros(4);
dW_x1_mat = zeros(4);
dW_x2_mat = zeros(4);
W_upper_mat = zeros(4);

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');
        dec_sol = clean(double(dec_List),1e-3);
        
        %         W_sol = replace(W,dec_List,dec_sol);
        %
        %         dW_x1 = jacobian(W(:),x1);
        %         dW_x1_sol = replace(reshape(dW_x1,n,n),dec_List,dec_sol);
        %
        %         dW_x2 = jacobian(W(:),x2);
        %         dW_x2_sol = replace(reshape(dW_x2,n,n),dec_List,dec_sol);
        
        W_sol = zeros(n,n,length(w_mono));
        for i = 1:length(w_mono)
            W_sol(:,:,i) = replace(W_list{i},dec_List,dec_sol);
        end
        
        dw_mono_x1 = jacobian(w_mono,x1);
        dw_mono_x2 = jacobian(w_mono,x2);
        
        W_upper_mat = replace(W_upper,dec_List,dec_sol);
        
        sdisplay(W_sol)
        pause;
        
        %% Create monomial functions
        for i = 1:length(w_mono)
            w_mono_i = sdisplay(w_mono(i));
            dw_mono_x1_i = sdisplay(dw_mono_x1(i));
            dw_mono_x2_i = sdisplay(dw_mono_x2(i));
            
            w_mono_i = strcat('@(x1)',w_mono_i{1});
            dw_mono_x1_i = strcat('@(x1)',dw_mono_x1_i{1});
            dw_mono_x2_i = strcat('@(x1,x2)',dw_mono_x2_i{1});
            
            eval(sprintf('w_mono_%d_fnc = str2func(w_mono_i);',i));
            eval(sprintf('dw_mono_x1_%d_fnc = str2func(dw_mono_x1_i);',i));
            eval(sprintf('dw_mono_x2_%d_fnc = str2func(dw_mono_x2_i);',i));
        end
        
        %% Create functions
        
        %         for i = 1:4
        %             for j = 1:4
        %                 W_ij = sdisplay(W_sol(i,j));
        %                 dW_x1_ij = sdisplay(dW_x1_sol(i,j));
        %                 dW_x2_ij = sdisplay(dW_x2_sol(i,j));
        %
        %                 W_ij = strcat('@(x1)',W_ij{1});
        %                 dW_x1_ij = strcat('@(x1)',dW_x1_ij{1});
        %                 dW_x2_ij = strcat('@(x1,x2)',dW_x2_ij{1});
        %
        %                 eval(sprintf('W_%d%d_fnc = str2func(W_ij);',i,j));
        %                 eval(sprintf('dW_x1_%d%d_fnc = str2func(dW_x1_ij);',i,j));
        %                 eval(sprintf('dW_x2_%d%d_fnc = str2func(dW_x2_ij);',i,j));
        %             end
        %         end
        
        %% Put together
        W_exec = 'W_mat = @(x)';
        dW_x1_exec = 'dW_x1_mat = @(x)';
        dW_x2_exec = 'dW_x2_mat = @(x)';
        
        for i = 1:length(w_mono)
            if i<length(w_mono)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*w_mono_%d_fnc(x(1)) +',i,i));
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*w_mono_%d_fnc(x(1));',i,i));
            end
        end
        for i = 2:length(w_mono)
            if i<length(w_mono)
                dW_x1_exec = strcat(dW_x1_exec,sprintf('W_sol(:,:,%d)*dw_mono_x1_%d_fnc(x(1)) +',i,i));
                dW_x2_exec = strcat(dW_x2_exec,sprintf('W_sol(:,:,%d)*dw_mono_x2_%d_fnc(x(1),x(2)) +',i,i));
            else
                dW_x1_exec = strcat(dW_x1_exec,sprintf('W_sol(:,:,%d)*dw_mono_x1_%d_fnc(x(1));',i,i));
                dW_x2_exec = strcat(dW_x2_exec,sprintf('W_sol(:,:,%d)*dw_mono_x2_%d_fnc(x(1),x(2));',i,i));
            end
        end
            
        %% Assemble
        
        %         W_exec = 'W_mat = @(x)[';
        %         dW_x1_exec = 'dW_x1_mat = @(x)[';
        %         dW_x2_exec = 'dW_x2_mat = @(x)[';
        %
        %         for i = 1:4
        %             for j = 1:4
        %                 if j<n
        %                     W_exec = strcat(W_exec,sprintf('W_%d%d_fnc(x(1)),',i,j));
        %                     dW_x1_exec = strcat(dW_x1_exec,sprintf('dW_x1_%d%d_fnc(x(1)),',i,j));
        %                     dW_x2_exec = strcat(dW_x2_exec,sprintf('dW_x2_%d%d_fnc(x(1),x(2)),',i,j));
        %                 elseif j==n
        %                     if i<n
        %                         W_exec = strcat(W_exec,sprintf('W_%d%d_fnc(x(1));',i,j));
        %                         dW_x1_exec = strcat(dW_x1_exec,sprintf('dW_x1_%d%d_fnc(x(1));',i,j));
        %                         dW_x2_exec = strcat(dW_x2_exec,sprintf('dW_x2_%d%d_fnc(x(1),x(2));',i,j));
        %                     else
        %                         W_exec = strcat(W_exec,sprintf('W_%d%d_fnc(x(1))];',i,j));
        %                         dW_x1_exec = strcat(dW_x1_exec,sprintf('dW_x1_%d%d_fnc(x(1))];',i,j));
        %                         dW_x2_exec = strcat(dW_x2_exec,sprintf('dW_x2_%d%d_fnc(x(1),x(2))];',i,j));
        %                     end
        %                 end
        %             end
        %         end
        
        %% Execute
        eval(W_exec); eval(dW_x1_exec); eval(dW_x2_exec);
        
    end
end
end