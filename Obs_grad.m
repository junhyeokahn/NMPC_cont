function J_grad_obs = Obs_grad(xu,n,m,N,obs)
%Dynamics, and terminal

global US_B;

J_grad_obs = zeros((N+1)*(n+m),1);

for i = 1:obs.n_obs
    o_pos = obs.pos(:,i);
    o_r = obs.r(i) + obs.infl; 
    decay = obs.decay(i);
    
    for k = 1:N+1
        J_grad_obs(1+(k-1)*n:2+(k-1)*n) = J_grad_obs(1+(k-1)*n:2+(k-1)*n) + ...
                                          US_B(k,i)*(-(decay*2/(o_r^2))*(xu(1+(k-1)*n:2+(k-1)*n)-o_pos));
    end
end

end