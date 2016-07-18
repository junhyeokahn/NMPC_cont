function conJ = NMPC_conJM(xu,Prob)
%Dynamics
n = Prob.user.n;
N = Prob.user.N;
m = Prob.user.m;

conJ = zeros(n*(N+1),(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tf)*Prob.user.D - ...
            df_all(Prob.user.df,xu(1:n*(N+1)),n,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

end

