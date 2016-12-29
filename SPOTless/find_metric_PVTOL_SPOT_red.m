function [solved,w_lower,w_upper] = ...
    find_metric_PVTOL_SPOT_red(n,g,p_lim,pd_lim,vy_lim,vz_lim,...
    condn,lambda,ccm_eps,return_metric)
%%


% W_scale = diag([0.02;0.0268;0.0367;0.0089;0.02;0.02]);
W_scale = diag([0.001;0.02;0.005;0.001;0.0367]);
% W_scale = zeros(n);

% sin_x = @(x) 0.05059*(x/(pi/6));
% cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

% sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
% cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

%states
x = msspoly('x',5);

%pos_def indeterminates
dfiv = msspoly('dfiv',5);
dtre = msspoly('dtre',3);

sin_p = sin_x(x(5));
cos_p = cos_x(x(5));

%dynamics f
f = [x(3)*cos_p - x(4)*sin_p;
    x(3)*sin_p + x(4)*cos_p;
    -g*sin_p;
    -g*cos_p;
    0];

%          Y Z  x(3)    x(4)    p
df =      [0,0,cos_p,-sin_p,-x(3)*sin_p-x(4)*cos_p;
           0,0,sin_p,cos_p,x(3)*cos_p-x(4)*sin_p;
           0,0,0,0,-g*cos_p;
           0,0,0,0, g*sin_p;
           zeros(1,5)];

B_perp = [eye(2), zeros(2,1);
          zeros(3,2), [1;0;-x(4)]];

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dfiv);
prog = prog.withIndeterminate(dtre);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);

%% Parametrize W (2)

w_poly = monomials([x(3), x(5)], 0:4);
W_list = cell(length(w_poly),1);
[prog, W_list{1}] = prog.newSym(5);
W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    [prog, W_list{i}] = prog.newSym(5);
    W = W + W_list{i}*w_poly(i);
end
    
dW_f = diff(W(:),x)*f;
dW_f = reshape(dW_f,5,5);

[prog, W_upper] = prog.newSym(5);

%% Killing field condns
% S = [0, 1;
%     -1, 0];
% 
% W_12 = W(1:2,3:4);
% W_22 = W(3:4,3:4);
% W_2c = W(3:4,5);
% 
% kill_12 = reshape(diff(W_12(:),[x(3);x(5)])*[x(4);1],2,2) - (W_12*S');
% kill_22 = reshape(diff(W_22(:),[x(3);x(5)])*[x(4);1],2,2) - (S*W_22 + W_22*S');
% kill_2c = diff(W_2c,[x(3);x(5)])*[x(4);1] - (S*W_2c);

% prog = prog.withPolyEqs([kill_12(:);
%                          kill_22(:);
%                          kill_2c]);

b2 = [zeros(2,1);x(4);-x(3);1];
db2_dx = diff(b2,x);
dW_b2 = diff(W(:),x)*b2;
dW_b2 = reshape(dW_b2,5,5);
kill = B_perp'*(dW_b2 - (db2_dx*W + W*db2_dx'))*B_perp;
prog = prog.withPolyEqs(kill(:));

%%

%Lagrange multipliers
box_lim = [ p_lim^2-x(5)^2;
            vy_lim^2-x(3)^2;
            vz_lim^2-x(4)^2];

l_order = 4;
l_def_states = [x(3);x(5); dfiv]';

[prog, Ll] = prog.newSOSPoly(monomials(l_def_states,0:l_order),2);
[prog, Lu] = prog.newSOSPoly(monomials(l_def_states,0:l_order),2);

l_ccm_states = [x(3);x(5);x(4);dtre]';
lc_order = 4;

[prog, Lc_1] = prog.newSOSPoly(monomials(l_ccm_states,0:lc_order+2),1);
[prog, Lc_2_3] = prog.newSOSPoly(monomials(l_ccm_states,0:lc_order),2);
Lc = [Lc_1; Lc_2_3];

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
% prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSOS((dfiv'*W*dfiv - w_lower*(dfiv'*dfiv)) - Ll'*box_lim(1:2));
prog = prog.withSOS(dfiv'*(W_upper - W)*dfiv - Lu'*box_lim(1:2));

%CCM condition
R_CCM = -B_perp'*(-dW_f + df*W + W*df' + 2*lambda*W)*B_perp;
prog = prog.withSOS((dtre'*R_CCM*dtre - ccm_eps*(dtre'*dtre)) - Lc'*box_lim);

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.coneVar; prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

SOS_soln = prog.minimize(trace(W_scale*W_upper) + (1e-3)*sum(a), @spot_mosek, options);

solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

%% Parse

w_lower = double(SOS_soln.eval(w_lower));
w_upper = double(SOS_soln.eval(w_upper));


if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-3);
        end
        
        dw_poly_p = diff(w_poly,x(5));
        dw_poly_vy = diff(w_poly,x(3));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-3);
        
%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_p_fnc = mss2fnc(dw_poly_p,x,randn(length(x),2));
        dw_poly_vy_fnc = mss2fnc(dw_poly_vy,x,randn(length(x),2));
        
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
        save('metric_PVTOL_red_vectorized.mat','W_eval','w_poly_fnc','dw_poly_p_fnc','dw_poly_vy_fnc','W_upper');
    end
end
end