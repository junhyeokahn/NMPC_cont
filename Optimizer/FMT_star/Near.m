function N = Near(V,q,EPS)

%N = {u in V\{q} : J(q,u) <= r} 
N = zeros(length(V),1);
count = 1;
for i = 1:length(V)
    if (V(i).idx~=q.idx) && (dist_3d(V(i).coord,q.coord)<= EPS.r)
        N(count) = V(i).idx;
        count = count+1;
    end
end
N = N(1:count-1);

end