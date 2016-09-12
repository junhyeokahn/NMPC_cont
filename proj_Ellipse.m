function proj_Ellipse(dim, A, bound,center, N, c)

n = size(A,1);

In = eye(n);
P = In(dim,:);
Ap = (P*(A\P'))\eye(2);
Ellipse_plot(Ap*(1/bound),center,N,c);

end