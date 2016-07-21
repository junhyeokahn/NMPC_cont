function aux_Prob = setup_opt_aux(m,eps_u)

Name = 'Aux Problem';
% aux_Prob = qpAssign(Q_bar,zeros(1,(K_e+1)*m),ones(K_e+1,(K_e+1)*m),...
%          -Inf*ones((K_e+1)*m,1),zeros((K_e+1)*m,1),[],[],...
%          zeros((K_e+1)*m,1),Name);

aux_Prob = qpAssign((1+eps_u)*eye(m),zeros(m,1),ones(1,m),...
         -Inf,0,[],[],...
         zeros(m,1),Name);
       


end