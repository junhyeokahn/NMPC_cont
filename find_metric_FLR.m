function [solved, w_lower, w_upper, W_mat, dW_x1_mat, dW_x2_mat] = find_metric_FLR(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
                    condn, lambda, ccm_eps,return_metric)
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
W = sdpvar(n);
for i = 1:n
    for j = 1:n
        if (j>=i)
            [W(i,j),cp_c] = polynomial(w_states,w_order);
            dec_List = [dec_List;cp_c];
        else
            W(i,j) = W(j,i);
        end
    end
end

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
W_bounds = [w_lower>=1, w_upper >= w_lower];

%Condition bound
W_cond = [w_upper-condn*w_lower <= 0];

%W pos def
p_low = (delta_4'*W*delta_4 - w_lower*(delta_4'*delta_4)) - Ll*box_lim(1:2);
p_high = (w_upper*(delta_4'*delta_4) - delta_4'*W*delta_4) - Lu*box_lim(1:2);

%CCM condition
R_CCM = -(-dW_f(1:3,1:3) + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W(1:3,1:3));
p_CCM = (delta_3'*R_CCM*delta_3 - ccm_eps*(delta_3'*delta_3)) - Lc*box_lim;

%Assemble
constr_List = [W_bounds; W_cond; sos(p_low); sos(p_high); sos(p_CCM);
    sos(Ll_1); sos(Ll_2);% sos(Ll_3); sos(Ll_4);
    sos(Lu_1); sos(Lu_2);% sos(Lu_3); sos(Lu_4);
    sos(Lc_1); sos(Lc_2); sos(Lc_3); sos(Lc_4)];% sos(Lc_5); sos(Lc_6)];

options = sdpsettings('solver','mosek','verbose',0);
SOS_soln = solvesos(constr_List,(1e-4)*norm(dec_List,1),options, dec_List);

solved = SOS_soln.problem;
w_lower = double(w_lower);
w_upper = double(w_upper);

W_mat = zeros(4);
dW_x1_mat = zeros(4);
dW_x2_mat = zeros(4);

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');
        dec_sol = clean(double(dec_List),1e-3);
        W_sol = replace(W,dec_List,dec_sol);
        
        dW_x1 = jacobian(W(:),x1);
        dW_x1_sol = replace(reshape(dW_x1,n,n),dec_List,dec_sol);
        
        dW_x2 = jacobian(W(:),x2);
        dW_x2_sol = replace(reshape(dW_x2,n,n),dec_List,dec_sol);
        
        sdisplay(W_sol)
        pause;
                
        %% Create functions

        W_cell = cell(4,4);
        dW_x1_cell = cell(4,4);
%         dW_x2_cell = cell(4,4);
        
        for i = 1:4
            for j = 1:4
                W_ij = sdisplay(W_sol(i,j));
                dW_x1_ij = sdisplay(dW_x1_sol(i,j));
%                 dW_x2_ij = sdisplay(dW_x2_sol(i,j));
                
                W_ij = strcat('@(x1)',W_ij{1});
                dW_x1_ij = strcat('@(x1)',dW_x1_ij{1});
%                 dW_x2_ij = strcat('@(x1,x2)',dW_x2_ij{1});
                
                W_cell{i,j} = str2func(W_ij);
                dW_x1_cell{i,j} = str2func(dW_x1_ij);
%                 dW_x2_cell{i,j} = str2func(dW_x2_ij);
            end
        end
        W_mat = @(x) [W_cell{1,1}(x(1)), W_cell{1,2}(x(1)), W_cell{1,3}(x(1)), W_cell{1,4}(x(1));
                      W_cell{2,1}(x(1)), W_cell{2,2}(x(1)), W_cell{2,3}(x(1)), W_cell{2,4}(x(1));
                      W_cell{3,1}(x(1)), W_cell{3,2}(x(1)), W_cell{3,3}(x(1)), W_cell{3,4}(x(1));
                      W_cell{4,1}(x(1)), W_cell{4,2}(x(1)), W_cell{4,3}(x(1)), W_cell{4,4}(x(1))];
        
        dW_x1_mat = @(x) [dW_x1_cell{1,1}(x(1)), dW_x1_cell{1,2}(x(1)), dW_x1_cell{1,3}(x(1)), dW_x1_cell{1,4}(x(1));
                          dW_x1_cell{2,1}(x(1)), dW_x1_cell{2,2}(x(1)), dW_x1_cell{2,3}(x(1)), dW_x1_cell{2,4}(x(1));
                          dW_x1_cell{3,1}(x(1)), dW_x1_cell{3,2}(x(1)), dW_x1_cell{3,3}(x(1)), dW_x1_cell{3,4}(x(1));
                          dW_x1_cell{4,1}(x(1)), dW_x1_cell{4,2}(x(1)), dW_x1_cell{4,3}(x(1)), dW_x1_cell{4,4}(x(1))];
        
%         dW_x2_mat = @(x) [dW_x2_cell{1,1}(x(1)), dW_x2_cell{1,2}(x(1)), dW_x2_cell{1,3}(x(1)), dW_x2_cell{1,4}(x(1));
%                           dW_x2_cell{2,1}(x(1)), dW_x2_cell{2,2}(x(1)), dW_x2_cell{2,3}(x(1)), dW_x2_cell{2,4}(x(1));
%                           dW_x2_cell{3,1}(x(1)), dW_x2_cell{3,2}(x(1)), dW_x2_cell{3,3}(x(1)), dW_x2_cell{3,4}(x(1));
%                           dW_x2_cell{4,1}(x(1)), dW_x2_cell{4,2}(x(1)), dW_x2_cell{4,3}(x(1)), dW_x2_cell{4,4}(x(1))];
    end
end
end