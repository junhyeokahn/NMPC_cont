function c = NMPC_conMN(xu,Prob) 
%Dynamics, and terminal

n = Prob.user.n;
N = Prob.user.N;

c = zeros(n*(N+1)+1,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D*xu(1:n*(N+1)) -...
    (NMPC_dyn(Prob.user.f,xu(1:n*(N+1)),n,N) + Prob.user.B_full*xu(n*(N+1)+1:end));

%% Terminal constraint
c(end) = xu(n*N+1:n*(N+1))'*Prob.user.P*(xu(n*N+1:n*(N+1)));

end

